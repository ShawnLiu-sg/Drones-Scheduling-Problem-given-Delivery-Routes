"""Independent time-indexed reference MILP and physical witness verification.

One dispatch binary is indexed by (team, stop, integer start). Covering every
tick of its half-open round trip enforces the capacity of the shared depot
fleet. This formulation uses neither Benders certificates nor schedule columns.
"""

from collections import Counter, defaultdict
from math import isfinite
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB


def _latest(route, parcel):
    if parcel["kind"] == "truck":
        return parcel["dest"]
    if parcel["kind"] == "sortie":
        return route["sorties"][parcel["dest"]]["l"]
    raise ValueError(f"Unknown parcel kind {parcel['kind']}")


def _clean(x):
    return int(round(x)) if abs(x - round(x)) <= 1e-7 else x


def validate(I, routes):
    """Validate the original, joint schedule; return sum of team return delays.

    Fleet concurrency is checked at actual half-open interval endpoints, not
    sampled on a grid. Integer starts are required only by the reference solver.
    """
    tol = 1e-7

    def require(condition, msg):
        if not condition:
            raise ValueError(msg)

    require(routes is not None and len(routes) == len(I["routes"]), "route count")
    fleet_changes = Counter()
    total = 0.0
    for k, (route, plan) in enumerate(zip(I["routes"], routes)):
        n = len(route["service"])
        sorties = route["sorties"]
        parcels = route["parcels"]
        prefix = f"route {k}: "
        A, D, L, U = (plan[s] for s in ("A", "D", "L", "U"))
        assign = plan["assign"]
        starts = {int(i): t for i, t in plan["starts"].items()}
        require(len(A) == len(D) == n and len(L) == len(U) == len(sorties),
                prefix + "event counts")
        require(len(assign) == len(parcels), prefix + "parcel count")
        times = [*A, *D, *L, *U, *starts.values()]
        require(all(isfinite(t) and -tol <= t <= I["H"] + tol for t in times),
                prefix + "event outside horizon")
        require(abs(A[0]) <= tol and abs(D[0]) <= tol, prefix + "depot departure")
        for i in range(1, n):
            require(abs(A[i] - (D[i - 1] + route["travel"][i - 1])) <= tol,
                    prefix + f"travel to stop {i}")
            require(D[i] + tol >= A[i] + route["service"][i],
                    prefix + f"service at stop {i}")
        for s, sortie in enumerate(sorties):
            i, j = sortie["l"], sortie["b"]
            require(A[i] - tol <= L[s] <= D[i] + tol, prefix + f"launch {s}")
            require(A[j] - tol <= U[s] <= D[j] + tol, prefix + f"recovery {s}")
            require(sortie["F"] - tol <= U[s] - L[s] <= sortie["E"] + tol,
                    prefix + f"flight duration {s}")
            require(s == 0 or L[s] + tol >= U[s - 1], prefix + "sortie overlap")
        require(set(assign) == set(starts), prefix + "occupied stops differ from parcel stops")
        load = Counter()
        for p, i in zip(parcels, assign):
            require(isinstance(i, int) and 0 < i < n - 1 and route["tau"][i] > 0,
                    prefix + "ineligible stop")
            require(i <= _latest(route, p), prefix + "parcel arrived too late on route")
            require(starts[i] + tol >= p["r"], prefix + "dispatch before release")
            h = starts[i] + route["tau"][i]
            if p["kind"] == "truck" and p["dest"] == i:
                require(D[i] + tol >= h + route["service"][i],
                        prefix + "same-stop parcel not available for service")
            if p["kind"] == "sortie":
                require(L[p["dest"]] + tol >= h, prefix + "flight departed before receipt")
            load[i] += p["w"]
        require(all(w <= I["Q"] + tol for w in load.values()), prefix + "payload")
        for i, t in starts.items():
            tau = route["tau"][i]
            require(A[i] - tol <= t + tau <= D[i] + tol, prefix + "handoff window")
            require(t + 2 * tau <= I["H"] + tol, prefix + "round trip outside horizon")
            fleet_changes[t] += 1
            fleet_changes[t + 2 * tau] -= 1
        require(D[-1] + tol >= route["base"], prefix + "completion below baseline")
        total += D[-1] - route["base"]
    active = 0
    for t, delta in sorted(fleet_changes.items()):
        active += delta
        require(active <= I["M"], f"shared fleet overload at {t}: {active}")
    return _clean(total)


def solve(I, time_limit=30, threads=1, *, fixed_assignments=None, progress=None):
    """Solve the finite integer-time instance; integer data admit an optimal grid schedule."""
    begin = perf_counter()
    H = I["H"]
    if not isinstance(H, int):
        raise ValueError("H must be expressed in integer ticks")
    warm = I.get("warm")
    warm_obj = validate(I, warm) if warm is not None else None
    if fixed_assignments is not None and warm is not None:
        if [list(w["assign"]) for w in warm] != [list(a) for a in fixed_assignments]:
            raise ValueError("Common initial schedule violates fixed policy")
    if progress is not None:
        progress(warm_obj,0.)
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", 0)
        env.start()
        with gp.Model("reference", env=env) as model:
            model.Params.Threads = threads
            model.Params.MIPGap = 0
            model.Params.MIPGapAbs = 0
            model.Params.Seed = 13
            model.Params.FeasibilityTol = 1e-8
            model.Params.IntFeasTol = 1e-8
            model.Params.OptimalityTol = 1e-8
            model.ModelSense = GRB.MINIMIZE
            occupancy = defaultdict(list)
            vars_by_route = []
            terminals = []
            for k, r in enumerate(I["routes"]):
                n, m = len(r["service"]), len(r["sorties"])
                parcels = r["parcels"]
                duration_ub = min(H, int(warm_obj + r["base"])) if warm is not None else H
                A, D, L, U = (model.addVars(count, ub=duration_ub, name=f"{sym}_{k}")
                              for sym, count in (("A", n), ("D", n), ("L", m), ("U", m)))
                model.addConstr(A[0] == 0)
                model.addConstr(D[0] == 0)
                for i in range(1, n):
                    model.addConstr(A[i] == D[i-1] + r["travel"][i-1])
                    model.addConstr(D[i] >= A[i] + r["service"][i])
                for s, so in enumerate(r["sorties"]):
                    i, j = so["l"], so["b"]
                    model.addConstr(L[s] >= A[i])
                    model.addConstr(L[s] <= D[i])
                    model.addConstr(U[s] >= A[j])
                    model.addConstr(U[s] <= D[j])
                    model.addConstr(U[s] - L[s] >= so["F"])
                    model.addConstr(U[s] - L[s] <= so["E"])
                    if s:
                        model.addConstr(L[s] >= U[s-1])
                eligible = {j: [i for i in range(1, n-1)
                                if r["tau"][i] > 0 and i <= _latest(r, p)]
                            for j, p in enumerate(parcels)}
                if fixed_assignments is not None:
                    assign=fixed_assignments[k]
                    if len(assign)!=len(parcels) or any(i not in eligible[j] for j,i in enumerate(assign)):
                        raise ValueError("Ineligible fixed handoff")
                    eligible={j:[i] for j,i in enumerate(assign)}
                x = model.addVars([(j, i) for j in eligible for i in eligible[j]],
                                  vtype=GRB.BINARY, name=f"parcel_{k}")
                for j, stops in eligible.items():
                    model.addConstr(gp.quicksum(x[j, i] for i in stops) == 1)
                y = {}
                for i in sorted({i for stops in eligible.values() for i in stops}):
                    tau = r["tau"][i]
                    choices = [j for j in eligible if i in eligible[j]]
                    earliest = min(parcels[j]["r"] for j in choices)
                    # H bounds final returns. Local team completion bounds the
                    # handoff, even when the depot return leg extends later.
                    latest = min(H - 2*tau, duration_ub - tau)
                    yi = model.addVars(range(earliest, latest + 1),
                                       vtype=GRB.BINARY, name=f"depart_{k}_{i}")
                    y[i] = yi
                    active = yi.sum()
                    start = gp.quicksum(t * v for t, v in yi.items())
                    handoff = start + tau * active
                    model.addConstr(active <= 1)
                    model.addConstr(active <= gp.quicksum(x[j, i] for j in choices))
                    model.addConstr(gp.quicksum(parcels[j]["w"] * x[j, i] for j in choices)
                                    <= I["Q"] * active)
                    model.addConstr(A[i] <= handoff + duration_ub * (1 - active))
                    model.addConstr(D[i] >= handoff)
                    for j in choices:
                        p = parcels[j]
                        model.addConstr(x[j, i] <= active)
                        model.addConstr(start >= p["r"] * x[j, i])
                        if p["kind"] == "truck" and p["dest"] == i:
                            model.addConstr(D[i] >= handoff + r["service"][i]*x[j, i])
                        if p["kind"] == "sortie":
                            model.addConstr(L[p["dest"]] >=
                                            handoff - duration_ub * (1 - x[j, i]))
                    for t, v in yi.items():
                        for u in range(t, t + 2*tau):
                            occupancy[u].append(v)
                vars_by_route.append((A, D, L, U, x, y, eligible))
                terminals.append(D[n-1] - r["base"])
                if warm is not None:
                    w = warm[k]
                    for vs, name in ((A, "A"), (D, "D"), (L, "L"), (U, "U")):
                        for j, v in vs.items():
                            v.Start = w[name][j]
                    for (j, i), v in x.items():
                        v.Start = int(w["assign"][j] == i)
                    starts = {int(i): t for i, t in w["starts"].items()}
                    for i, yi in y.items():
                        for t, v in yi.items():
                            v.Start = int(starts.get(i) == t)
            for t, intervals in occupancy.items():
                model.addConstr(gp.quicksum(intervals) <= I["M"], name=f"fleet_{t}")
            model.setObjective(gp.quicksum(terminals))
            model.Params.TimeLimit = max(0.0, time_limit - (perf_counter() - begin))
            def callback(m,where):
                if where==GRB.Callback.MIP:
                    progress(m.cbGet(GRB.Callback.MIP_OBJBST),m.cbGet(GRB.Callback.MIP_OBJBND))
                elif where==GRB.Callback.MIPSOL:
                    progress(m.cbGet(GRB.Callback.MIPSOL_OBJ),m.cbGet(GRB.Callback.MIPSOL_OBJBND))
            if progress is None:
                model.optimize()
            else:
                model.optimize(callback)
            status = "optimal" if model.Status == GRB.OPTIMAL else (
                "infeasible" if model.Status == GRB.INFEASIBLE else "limit")
            routes = None
            obj = None
            if model.SolCount:
                routes = []
                for A,D,L,U,x,y,eligible in vars_by_route:
                    z = {s: [_clean(v.X) for v in vs.values()]
                         for s, vs in (("A", A), ("D", D), ("L", L), ("U", U))}
                    z["assign"] = [next(i for i in eligible[j] if x[j, i].X > .5)
                                   for j in eligible]
                    z["starts"] = {i: t for i, yi in y.items()
                                   for t, v in yi.items() if v.X > .5}
                    routes.append(z)
                obj = validate(I, routes)
                if abs(obj - model.ObjVal) > 1e-6:
                    raise RuntimeError("objective disagrees with physical witness")
            if warm is not None and (obj is None or warm_obj < obj):
                if status == "infeasible":
                    raise RuntimeError("reference declared a warm-feasible instance infeasible")
                obj, routes = warm_obj, warm
            lb = None
            if status == "optimal":
                lb = obj
            elif status != "infeasible":
                lb = max(0.0, model.ObjBound) if abs(model.ObjBound) < GRB.INFINITY else 0.0
            if progress is not None:
                progress(obj,lb)
            return dict(status=status, objective=obj, bound=lb, routes=routes,
                        stats=dict(seconds=perf_counter() - begin, nodes=model.NodeCount,
                                   variables=model.NumVars, binaries=model.NumBinVars,
                                   constraints=model.NumConstrs, solver_status=model.Status))
