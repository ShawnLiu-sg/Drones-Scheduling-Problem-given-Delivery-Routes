"""Exact branch-and-price on complete, timed single-team columns.

At every node we solve the master by column generation. Pricing uses a separate
MILP for each team. Branching fixes a stop/start-time binary in that team's
pricing MILP and filters inherited columns. The common integer unit is exact
for integer difference lags: fix any feasible set of drone chains and take the
earliest schedule, whose event times are integer path lengths.
"""
import heapq
import math
import time

import gurobipy as gp
from gurobipy import GRB

from .model import eligible, realize


def _model(env, threads):
    model = gp.Model(env=env)
    model.Params.OutputFlag = 0
    model.Params.Threads = threads
    model.Params.Seed = 13
    model.Params.MIPGap = 0
    model.Params.MIPGapAbs = 0
    model.Params.FeasibilityTol = 1e-8
    model.Params.OptimalityTol = 1e-8
    model.Params.IntFeasTol = 1e-8
    return model


class _Price:
    def __init__(self, route, capacity, horizon, env, threads, flight_horizon=None):
        self.flight_horizon = horizon+max(route['tau']) if flight_horizon is None else flight_horizon
        self.rt = route
        self.m = m = _model(env, threads)
        self.horizon = horizon
        self.base = route["base"]
        n = len(route["service"])
        a = m.addVars(n, ub=horizon, name="arrival")
        d = m.addVars(n, ub=horizon, name="departure")
        a[0].UB = d[0].UB = 0
        for i in range(1, n):
            m.addConstr(a[i] == d[i - 1] + route["travel"][i - 1])
            m.addConstr(d[i] >= a[i] + route["service"][i])
        legs = route["sorties"]
        launch = m.addVars(len(legs), ub=horizon, name="launch")
        recovery = m.addVars(len(legs), ub=horizon, name="recovery")
        for s, leg in enumerate(legs):
            for event, i in ((launch[s], leg["l"]), (recovery[s], leg["b"])):
                m.addConstr(event >= a[i])
                m.addConstr(event <= d[i])
            m.addConstr(recovery[s] >= launch[s] + leg["F"])
            m.addConstr(launch[s] >= recovery[s] - leg["E"])
            if s:
                m.addConstr(launch[s] >= recovery[s - 1])
        ps = route["parcels"]
        self.allowed = allowed = [eligible(route, p) for p in ps]
        self.x = x = m.addVars([(j,i) for j in range(len(ps)) for i in allowed[j]],
                               vtype=GRB.BINARY, name="parcel")
        for j, opts in enumerate(allowed):
            m.addConstr(gp.quicksum(x[j,i] for i in opts) == 1)
        # Tighten the time index using mandatory truck travel/service only.
        amin = [0] * n
        suffix = [0] * n
        for i in range(1, n):
            amin[i] = amin[i-1] + route["service"][i-1] + route["travel"][i-1]
        for i in reversed(range(n-1)):
            suffix[i] = suffix[i+1] + route["travel"][i] + route["service"][i+1]
        self.y = {}
        for i, tau in enumerate(route["tau"]):
            jobs = [j for j, opts in enumerate(allowed) if i in opts]
            if not jobs:
                continue
            z = m.addVar(vtype=GRB.BINARY, name=f"used[{i}]")
            earliest = max(0, amin[i]-tau, min(ps[j]["r"] for j in jobs))
            latest = min(horizon - suffix[i] - tau, self.flight_horizon-2*tau)
            ys = {t: m.addVar(vtype=GRB.BINARY, name=f"start[{i},{t}]")
                  for t in range(earliest,latest+1)}
            self.y.update({(i,t): v for t,v in ys.items()})
            m.addConstr(gp.quicksum(ys.values()) == z)
            m.addConstr(gp.quicksum(x[j,i] for j in jobs) >= z)
            m.addConstr(gp.quicksum(ps[j]["w"]*x[j,i] for j in jobs) <= capacity*z)
            r = m.addVar(ub=horizon, name=f"release[{i}]")
            m.addConstr(r == gp.quicksum(t*v for t,v in ys.items()))
            m.addGenConstrIndicator(z,1,r+tau >= a[i])
            m.addGenConstrIndicator(z,1,r+tau <= d[i])
            for j in jobs:
                m.addConstr(x[j,i] <= z)
                m.addConstr(r >= ps[j]["r"]*x[j,i])
                if ps[j]["kind"] == "sortie":
                    m.addGenConstrIndicator(x[j,i],1,r+tau <= launch[ps[j]["dest"]])
                elif ps[j]["dest"] == i:
                    m.addGenConstrIndicator(x[j,i],1,r+tau+route["service"][i] <= d[i])
        self.end = d[n-1]
        self.end.LB = self.base
        m.update()
        self.fixed = {}

    def price(self, dual, penalty_prefix, phase, fixed, deadline, cut_costs=None):
        m = self.m
        for k in self.fixed.keys() - fixed.keys():
            self.y[k].LB, self.y[k].UB = 0, 1
        for k, b in fixed.items():
            if self.fixed.get(k) != b:
                self.y[k].LB, self.y[k].UB = b, b
        self.fixed = fixed
        # The same time atoms are charged for their full half-open occupancy.
        objective = ((self.end-self.base if phase == 2 else 0) - dual
                     + gp.quicksum((penalty_prefix[t+2*self.rt["tau"][i]]-penalty_prefix[t])*v
                                   for (i,t),v in self.y.items()))
        if cut_costs:
            objective += gp.quicksum(cost*self.y[atom] for atom,cost in cut_costs.items())
        m.setObjective(objective)
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            return None, -math.inf
        m.Params.TimeLimit = remaining
        m.optimize()
        if m.Status == GRB.INFEASIBLE:
            return None, math.inf
        bound = m.ObjBound if m.Status in (GRB.OPTIMAL,GRB.TIME_LIMIT,GRB.INTERRUPTED) else -math.inf
        if not m.SolCount:
            return None, bound
        assign = [next(i for i in opts if self.x[j,i].X > .5) for j,opts in enumerate(self.allowed)]
        starts = {i:t for (i,t),v in self.y.items() if v.X > .5}
        witness = realize(self.rt,assign,starts)
        if witness is None or witness["D"][-1] > self.horizon:
            raise ArithmeticError("pricing MILP returned an inconsistent timetable")
        rc = ((witness["D"][-1]-self.base) if phase == 2 else 0) - dual
        rc += sum(penalty_prefix[t+2*self.rt["tau"][i]]-penalty_prefix[t] for i,t in starts.items())
        if cut_costs:
            rc += sum(cut_costs.get(atom,0.) for atom in starts.items())
        return (witness, rc), bound


def solve(I, time_limit=30, threads=1, *, strengthening=None, route_price=None, progress=None):
    start = time.perf_counter()
    deadline = start + time_limit
    routes = I["routes"]
    kmax = len(routes)
    stats = dict(nodes=0,pricing_calls=0,master_solves=0,columns=0,branches=0,
                 phase1_pricing=0,root_bound=None,pricing_seconds=0.0,
                 restricted_mip_calls=0,restricted_mip_improvements=0,restricted_mip_seconds=0.0)
    incumbent = I.get("warm")
    ub = sum(w["D"][-1]-r["base"] for w,r in zip(incumbent,routes)) if incumbent is not None else math.inf
    stats.update(root_dual_lb=0.0, root_lp_exact=None)
    stats.update(overload_columns=0, generated_overload=0, node_certificates=[])
    def tick_global(lb):
        if progress is not None:
            progress(ub,lb)
    tick_global(0)
    def result(status, lb):
        tick_global(lb)
        if status=="limit" and math.isfinite(ub) and lb==ub:
            status="optimal"
        if strengthening is not None:
            stats.update(strengthening.stats)
        stats["wall_seconds"] = time.perf_counter()-start
        return dict(status=status, objective=ub if math.isfinite(ub) else None,
                    bound=lb if math.isfinite(lb) else None, routes=incumbent, stats=stats)
    if ub == 0:
        stats["root_lp_exact"] = 0.0
        stats["root_bound"] = 0
        return result("optimal", 0)
    # A witness with sum delay U permits D_end[k] <= base[k]+U for every
    # solution improving or tying that witness, since all delays are nonnegative.
    caps = [int(min(I["H"],r["base"]+ub)) for r in routes]
    slots = max(caps)+max(max(r["tau"]) for r in routes)
    allcols = [{} for _ in routes]
    def add(k,w):
        w = dict(w,starts={int(i):int(t) for i,t in w["starts"].items()})
        key = frozenset(w["starts"].items())
        delay = w["D"][-1]-routes[k]["base"]
        old = allcols[k].get(key)
        if old is not None and old[0] <= delay:
            return False
        occ = {}
        for i,t in key:
            for s in range(t,t+2*routes[k]["tau"][i]):
                occ[s] = occ.get(s,0)+1
        allcols[k][key] = (delay,occ,w)
        if old is None and any(value>I["M"] for value in occ.values()):
            stats["overload_columns"] += 1
        stats["columns"] += int(old is None)
        return True
    if incumbent is not None:
        for k,w in enumerate(incumbent):
            add(k,w)
    env = gp.Env(empty=True)
    env.setParam("OutputFlag",0)
    env.start()
    pricers = []
    last_heur_size = 0

    def restricted_incumbent():
        """Search the true columns from all nodes. Results only update the UB."""
        nonlocal incumbent, ub, last_heur_size
        total = sum(map(len,allcols))
        tick = time.perf_counter()
        if (total-last_heur_size < kmax or any(not pool for pool in allcols)
                or deadline-tick < 0.02):
            return
        last_heur_size = total
        stats["restricted_mip_calls"] += 1
        with _model(env,threads) as h:
            eq = [h.addConstr(gp.LinExpr()==1) for _ in routes]
            caprows = [h.addConstr(gp.LinExpr()<=I["M"]) for _ in range(slots)]
            vv = []
            for k,pool in enumerate(allcols):
                win = (frozenset((int(i),int(t)) for i,t in incumbent[k]["starts"].items())
                       if incumbent is not None else None)
                for key,(delay,occ,w) in pool.items():
                    gc = gp.Column([1]+list(occ.values()),[eq[k]]+[caprows[t] for t in occ])
                    v = h.addVar(obj=delay,vtype=GRB.BINARY,column=gc)
                    if win is not None:
                        v.Start = int(key==win)
                    vv.append((k,v,w))
            rem = deadline-time.perf_counter()
            if rem > 0:
                h.Params.TimeLimit = min(rem,0.08)
                h.optimize()
                if h.SolCount and h.ObjVal < ub-1e-6:
                    wits = [None]*kmax
                    for k,v,w in vv:
                        if v.X > .5:
                            wits[k] = w
                    ub = sum(w["D"][-1]-r["base"] for w,r in zip(wits,routes))
                    incumbent = wits
                    stats["restricted_mip_improvements"] += 1
        stats["restricted_mip_seconds"] += time.perf_counter()-tick

    try:
        for r,cap in zip(routes,caps):
            if time.perf_counter() >= deadline:
                return result("limit", 0)
            pricers.append(_Price(r,I["Q"],cap,env,threads,I['H']))
        if strengthening is not None:
            strengthening.prepare(I,pricers,slots)
        queue = [(0,0,())]
        ticket = 0
        active_lb = 0
        while queue:
            active_lb,_,branches = heapq.heappop(queue)
            if active_lb >= ub:
                continue
            if time.perf_counter() >= deadline:
                return result("limit", min(ub,active_lb))
            tick_global(min(active_lb,ub))
            stats["nodes"] += 1
            fixes = [{} for _ in routes]
            for k,i,t,b in branches:
                fixes[k][i,t] = b
            # If feasible columns already cover all teams without overload,
            # phase I immediately gives zero and no pricing is needed in it.
            phase = 1
            while time.perf_counter() < deadline:
                if phase == 2 and stats["master_solves"] % 3 == 0:
                    restricted_incumbent()
                    tick_global(min(ub,active_lb,queue[0][0] if queue else math.inf))
                    if active_lb >= ub:
                        phase = 0
                        break
                with _model(env,threads) as m:
                    m.Params.Method = 1
                    conv = [m.addConstr(gp.LinExpr() == 1) for _ in routes]
                    caprows = [m.addConstr(gp.LinExpr() <= I["M"]) for _ in range(slots)]
                    cuts = strengthening.cuts if strengthening is not None else ()
                    cutrows = [m.addConstr(gp.LinExpr() <= 1) for _ in cuts]
                    cols = {}
                    for k,known in enumerate(allcols):
                        for key,(delay,occ,w) in known.items():
                            if any((atom in key) != flag for atom,flag in fixes[k].items()):
                                continue
                            gc = gp.Column([1]+list(occ.values()),[conv[k]]+[caprows[t] for t in occ])
                            for cut,row in zip(cuts,cutrows):
                                if cut.hits(k,key):
                                    gc.addTerms(1.,row)
                            cols[k,key] = m.addVar(obj=(delay if phase==2 else 0),column=gc)
                        if phase == 1:
                            m.addVar(obj=1,column=gp.Column([1],[conv[k]]))
                    remaining = deadline-time.perf_counter()
                    if remaining <= 0:
                        break
                    m.Params.TimeLimit = remaining
                    m.optimize()
                    stats["master_solves"] += 1
                    if m.Status != GRB.OPTIMAL:
                        break
                    value = m.ObjVal
                    if phase == 1 and value <= 1e-8:
                        phase = 2
                        continue
                    if phase == 2 and not branches and strengthening is not None and strengthening.separates:
                        weights = [(k,key,v.X) for (k,key),v in cols.items() if v.X > 1e-8]
                        if strengthening.separate(weights, deadline=deadline):
                            continue
                    # Capacity duals are nonpositive for this minimization LP.
                    # Clipping numerical positives gives admissible multipliers.
                    penalty = [max(0.,-r.Pi) for r in caprows]
                    duals = [r.Pi for r in conv]
                    dualval = sum(duals)-I["M"]*sum(penalty)
                    atom_costs = [None]*kmax
                    if cutrows:
                        cut_prices = [max(0.,-r.Pi) for r in cutrows]
                        dualval -= sum(cut_prices)
                        atom_costs = strengthening.price_terms(cut_prices)
                    prefix = [0.]
                    for p in penalty:
                        prefix.append(prefix[-1]+p)
                    added = False
                    lbs = []
                    for k,p in enumerate(pricers):
                        tick = time.perf_counter()
                        if route_price is None:
                            cand, lb = p.price(duals[k],prefix,phase,fixes[k],deadline,atom_costs[k])
                        else:
                            cand, lb = route_price(k,p,duals[k],prefix,phase,fixes[k],deadline,atom_costs[k])
                        stats["pricing_seconds"] += time.perf_counter()-tick
                        stats["pricing_calls"] += 1
                        stats["phase1_pricing"] += int(phase==1)
                        if phase == 1:
                            lb = min(lb,1-duals[k])  # artificial column
                        lbs.append(lb)
                        if cand and cand[1] < -1e-7:
                            n_before=stats["overload_columns"]
                            added = add(k,cand[0]) or added
                            stats["generated_overload"] += stats["overload_columns"]-n_before
                        if time.perf_counter() >= deadline:
                            break
                    if len(lbs) < kmax or any(x == -math.inf for x in lbs):
                        break
                    # Integer-valued objectives allow ceiling a valid dual LB.
                    dual_lb = dualval + sum(min(0.,b) for b in lbs)
                    if phase == 1:
                        if dual_lb > 1e-7:
                            phase = 0  # certified infeasible node
                            break
                        if added:
                            continue
                        break
                    active_lb = max(active_lb, math.ceil(dual_lb-1e-6))
                    tick_global(min(ub,active_lb,queue[0][0] if queue else math.inf))
                    capture=I.get('_capture_nodes',0)
                    if capture:
                        limit=8 if capture is True else int(capture)
                        rows=stats['node_certificates']
                        bucket=(bool(branches),bool(cutrows and any(cut_prices)))
                        if sum(tuple(v['bucket'])==bucket for v in rows)<limit:
                            rows.append(dict(bucket=bucket,branches=branches,duals=duals,prefix=prefix,
                                cut_prices=list(cut_prices) if cutrows else [],
                                cuts=[{k:sorted(atoms) for k,atoms in c.atoms.items()} for c in cuts],
                                pricing_lbs=[v if math.isfinite(v) else ('+inf' if v>0 else '-inf') for v in lbs],
                                dual_lb=dual_lb,integer_lb=active_lb,team_caps=caps,
                                restricted_value=value))
                    if not branches:
                        stats["root_bound"] = active_lb
                        stats["root_dual_lb"] = max(stats["root_dual_lb"],dual_lb)
                        if abs(value-dual_lb) <= 1e-6:
                            stats["root_lp_exact"] = value
                    if active_lb >= ub:
                        phase = 0
                        break
                    if added:
                        continue
                    if any(b < -1e-6 for b in lbs):
                        break
                    weights = {(k,key):v.X for (k,key),v in cols.items() if v.X > 1e-8}
                    masses = {}
                    for (k,key),v in weights.items():
                        for i,t in key:
                            atom = (k,i,t)
                            masses[atom] = masses.get(atom,0)+v
                    frac = [(abs(v-.5),atom) for atom,v in masses.items() if 1e-7 < v < 1-1e-7]
                    if frac:
                        _,atom = min(frac)
                        for b in (0,1):
                            ticket += 1
                            heapq.heappush(queue,(active_lb,ticket,branches+(atom+(b,),)))
                        stats["branches"] += 1
                    else:
                        wits = []
                        for k in range(kmax):
                            key = max(((key,v) for (kk,key),v in weights.items() if kk==k),
                                      key=lambda kv:kv[1])[0]
                            wits.append(allcols[k][key][2])
                        obj = sum(w["D"][-1]-r["base"] for w,r in zip(wits,routes))
                        if active_lb < obj:
                            raise ArithmeticError("integer master primal and certified bound disagree")
                        if obj < ub:
                            ub,incumbent = obj,wits
                            tick_global(min(ub,active_lb,queue[0][0] if queue else math.inf))
                    phase = 0
                    break
            if phase:
                return result("limit", min(ub,active_lb,queue[0][0] if queue else math.inf))
        return result("optimal" if math.isfinite(ub) else "infeasible", ub)
    finally:
        for p in pricers:
            p.m.dispose()
        env.dispose()
