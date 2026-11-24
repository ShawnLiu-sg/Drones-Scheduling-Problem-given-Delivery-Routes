"""Exact Liu-style depot flow for fixed truck routes and delivery sorties.

Every late parcel must be handed over and the cost is total return delay.
Travel, service, release and flight times are on the same integer grid.
Replenishment handoff is instantaneous, with the truck at the handoff stop.
"""

from collections import defaultdict
from itertools import pairwise
from math import ceil, isfinite
from pathlib import Path
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB

from scheduling.study import model as schedules
from scheduling.study.reference import _latest, validate

TOL = 1e-6


def _road_tail(route):
    """Minimum truck travel and service after departure from each stop."""
    tail = [0] * len(route["service"])
    for i in range(len(tail) - 2, -1, -1):
        tail[i] = tail[i + 1] + route["travel"][i] + route["service"][i + 1]
    return tail


def _team_lower_bound(route, baseline, tail, latest, capacity):
    """Bound return delay by the best possible handoff for each parcel alone."""
    earliest_return = baseline["D"][-1]
    for parcel, last in zip(route["parcels"], latest):
        if parcel["w"] > capacity:
            return None
        sortie = route["sorties"][parcel["dest"]] if parcel["kind"] == "sortie" else None
        candidates = []
        for i in range(1, last + 1):
            tau = route["tau"][i]
            if tau <= 0:
                continue
            receipt = max(baseline["A"][i], parcel["r"] + tau)
            remaining = tail[i]
            if sortie is None and i == parcel["dest"]:
                remaining += route["service"][i]
            bound = receipt + remaining
            if sortie is not None:
                bound = max(bound, receipt + sortie["F"] + tail[sortie["b"]])
            candidates.append(bound)
        if not candidates:
            return None
        earliest_return = max(earliest_return, min(candidates))
    return earliest_return - route["base"]


def _departure_groups(parcels, eligible, first, last):
    """Successive intervals on which exactly the same parcels can depart."""
    cuts = sorted({first, last + 1} | {
        parcels[p]["r"] for p in eligible if first < parcels[p]["r"] <= last
    })
    for start, end in pairwise(cuts):
        yield start, end, [p for p in eligible if parcels[p]["r"] <= start]


def _add_stop(model, parcels, eligible, stop, earliest_arrival, latest_departure,
              tau, horizon, capacity, outgoing, incoming, warm, stats):
    first = max(0, earliest_arrival - tau, min(parcels[p]["r"] for p in eligible))
    last = min(latest_departure - tau, horizon - 2 * tau)
    if first > last:
        return None

    spans = []
    group_vars = []
    trip_vars = []
    receipt_times = []
    loads_here = defaultdict(list)
    warm_start = warm["starts"].get(stop) if warm is not None else None

    for start, end, allowed in _departure_groups(parcels, eligible, first, last):
        flights = list(model.addVars(end - start, vtype=GRB.BINARY).values())
        group = model.addVar(ub=1)
        model.addLConstr(gp.quicksum(flights) == group)
        # Only one handoff state at this stop can be selected. Within a group,
        # every departure admits the same parcels, so they share load variables.
        loads = model.addVars(allowed, vtype=GRB.BINARY)
        model.addLConstr(
            gp.LinExpr([parcels[p]["w"] for p in allowed],
                       [loads[p] for p in allowed]) <= capacity * group
        )
        model.addLConstr(loads.sum() >= group)
        for p in allowed:
            model.addLConstr(loads[p] <= group)
            loads_here[p].append(loads[p])
            if warm is not None:
                loads[p].Start = int(
                    warm_start is not None and start <= warm_start < end
                    and warm["assign"][p] == stop
                )

        for dispatch, flight in enumerate(flights, start):
            outgoing[dispatch].append(flight)
            incoming[dispatch + 2 * tau].append(flight)
            if warm is not None:
                flight.Start = int(dispatch == warm_start)
        spans.append((start, flights))
        group_vars.append(group)
        trip_vars.extend(flights)
        receipt_times.extend(range(start + tau, end + tau))
        stats["handoff_states"] += len(flights)
        stats["departure_groups"] += 1
        stats["parcel_variables"] += len(allowed)

    selected = gp.quicksum(group_vars)
    handoff = gp.LinExpr(receipt_times, trip_vars)
    model.addLConstr(selected <= 1)
    return spans, selected, handoff, loads_here


def _delivery_clock(model, route, baseline, departure_ub, warm):
    """Truck stop events and the onboard drone's fixed sortie sequence."""
    service = route["service"]
    sorties = route["sorties"]
    n, m = len(service), len(sorties)
    a = model.addVars(n, lb=baseline["A"], ub=[u - s for u, s in zip(departure_ub, service)])
    d = model.addVars(n, lb=baseline["D"], ub=departure_ub)
    l = model.addVars(m, lb=baseline["L"], ub=[
        min(departure_ub[f["l"]], departure_ub[f["b"]] - f["F"]) for f in sorties
    ])
    u = model.addVars(m, lb=baseline["U"], ub=[departure_ub[f["b"]] for f in sorties])
    model.addLConstr(a[0] == 0)
    model.addLConstr(d[0] == 0)
    for i in range(1, n):
        model.addLConstr(a[i] == d[i - 1] + route["travel"][i - 1])
        model.addLConstr(d[i] >= a[i] + service[i])
    for j, f in enumerate(sorties):
        model.addLConstr(l[j] >= a[f["l"]])
        model.addLConstr(l[j] <= d[f["l"]])
        model.addLConstr(u[j] >= a[f["b"]])
        model.addLConstr(u[j] <= d[f["b"]])
        model.addLConstr(u[j] - l[j] >= f["F"])
        model.addLConstr(u[j] - l[j] <= f["E"])
        if j:
            model.addLConstr(l[j] >= u[j - 1])
    if warm is not None:
        for name, variables in (("A", a), ("D", d), ("L", l), ("U", u)):
            for j, v in variables.items():
                v.Start = warm[name][j]
    return a, d, l


def _add_team(model, route, baseline, tail, latest, lower, upper, instance,
              outgoing, incoming, warm, deadline, stats):
    parcels = route["parcels"]
    choices = [[] for _ in parcels]
    stops = {}
    endcap = route["base"] + upper
    departure_ub = [endcap - remaining for remaining in tail]
    # With zero service and no sorties, equality of handoff and departure
    # gives the smaller cumulative-delay model used by the original comparison.
    use_events = bool(route["sorties"] or any(route["service"]))
    if use_events:
        a, d, l = _delivery_clock(model, route, baseline, departure_ub, warm)
    previous = None
    for stop in range(1, max(latest, default=0) + 1):
        if perf_counter() >= deadline:
            return None
        tau = route["tau"][stop]
        if tau <= 0:
            continue
        eligible = [p for p, last in enumerate(latest) if stop <= last]
        stop_cap = departure_ub[stop]
        built = _add_stop(
            model, parcels, eligible, stop, baseline["A"][stop], stop_cap, tau,
            instance["H"], instance["Q"], outgoing, incoming, warm, stats,
        )
        if built is None:
            continue
        spans, selected, handoff, loads_here = built
        stops[stop] = spans
        for p, loads in loads_here.items():
            choices[p].extend((stop, v) for v in loads)

        if use_events:
            model.addLConstr(a[stop] <= handoff + stop_cap * (1 - selected))
            model.addLConstr(d[stop] >= handoff)
            for p, loads in loads_here.items():
                parcel = parcels[p]
                received = gp.quicksum(loads)
                if parcel["kind"] == "sortie":
                    model.addLConstr(l[parcel["dest"]] >= handoff - stop_cap * (1 - received))
                elif parcel["dest"] == stop:
                    model.addLConstr(d[stop] >= handoff + route["service"][stop] * received)
        else:
            delay = handoff - baseline["A"][stop] * selected
            cumulative = model.addVar(ub=upper)
            model.addLConstr(delay <= cumulative)
            model.addLConstr(cumulative <= delay + upper * (1 - selected))
            if previous is not None:
                model.addLConstr(cumulative >= previous)
            previous = cumulative

    for options in choices:
        model.addLConstr(gp.quicksum(load for _, load in options) == 1)
    terminal = model.addVar(lb=lower, ub=upper, obj=1)
    if use_events:
        model.addLConstr(terminal == d[len(tail) - 1] - route["base"])
    elif previous is not None:
        model.addLConstr(terminal >= previous)
    return stops, choices


def _add_fleet_flow(model, outgoing, incoming, horizon, fleet, deadline, stats):
    """One unit of flow follows each drone through flights and depot waiting."""
    times = sorted({0, horizon, *outgoing, *incoming})
    wait = list(model.addVars(len(times) - 1, ub=fleet).values())
    stats["depot_time_nodes"] = len(times)
    stats["wait_arcs"] = len(wait)
    for j, time in enumerate(times):
        if perf_counter() >= deadline:
            return False
        balance = gp.quicksum(outgoing[time]) - gp.quicksum(incoming[time])
        if j < len(wait):
            balance += wait[j]
        if j > 0:
            balance -= wait[j - 1]
        supply = fleet * (int(time == 0) - int(time == horizon))
        model.addLConstr(balance == supply)
    return True


def _witness(instance, model, teams):
    assignments, starts = [], []
    for stops, choices in teams:
        assignments.append([
            next(stop for stop, parcel in options if parcel.X > 0.5)
            for options in choices
        ])
        starts.append({
            stop: first + j
            for stop, spans in stops.items()
            for first, flights in spans
            for j, used in enumerate(model.getAttr("X", flights))
            if used > 0.5
        })
    return schedules.timetable(instance["routes"], assignments, fixed=starts)


def _finite_bound(value):
    return isfinite(value) and abs(value) < GRB.INFINITY


def solve(instance, time_limit=120, threads=1, *, progress=None):
    """Return a validated schedule and a lower bound in the Chapter 5 format.

    The wall budget includes model generation. The optional progress function
    receives incumbent and bound updates in the same form as reference.solve.
    An optional grid-based initial schedule is read from instance["warm"].
    """
    begun = perf_counter()
    deadline = begun + time_limit
    routes = instance["routes"]

    stats = dict(handoff_states=0, departure_groups=0, parcel_variables=0,
                 depot_time_nodes=0, wait_arcs=0, build_seconds=0.0,
                 optimizer_seconds=0.0)
    initial = instance.get("warm")
    best, objective = None, None

    def result(status, bound, origin="algorithm"):
        if progress is not None:
            progress(objective, bound)
        return dict(status=status, objective=objective, bound=bound, routes=best,
                    origin=origin, stats=dict(stats, seconds=perf_counter() - begun))

    if "H" not in instance:
        # Archived inputs can omit the horizon when the parcel relaxation
        # already proves infeasibility. Recompute that certificate.
        from scheduling.relaxation import lower_bound as parcel_bound

        if not isfinite(parcel_bound(routes)):
            return result("infeasible", None, "common")
        raise ValueError("H is required when the parcel relaxation is feasible")
    horizon = instance["H"]
    if not isinstance(horizon, int) or horizon < 0:
        raise ValueError("H must be a nonnegative integer number of ticks")

    baselines = schedules.timetable(routes, [[] for _ in routes])
    if baselines is None:
        return result("infeasible", None, "common")
    if any(b["D"][-1] != r["base"] for b, r in zip(baselines, routes)):
        raise ValueError("Each base must be the earliest return on the fixed route without resupply")
    tails = [_road_tail(route) for route in routes]
    latest = [[min(_latest(r, p), len(r["tau"]) - 2) for p in r["parcels"]] for r in routes]
    lower = [_team_lower_bound(r, b, tail, last, instance["Q"])
             for r, b, tail, last in zip(routes, baselines, tails, latest)]
    if any(lb is None or r["base"] + lb > horizon for r, lb in zip(routes, lower)):
        return result("infeasible", None, "common")

    common_bound = sum(lower)
    if initial is not None:
        validate(instance, initial)
        starts = [{int(i): t for i, t in p["starts"].items()} for p in initial]
        if any(t != int(t) for start in starts for t in start.values()):
            raise ValueError("The initial drone schedule must use integer departure times")
        starts = [{i: int(t) for i, t in start.items()} for start in starts]
        # Remove discretionary truck waiting while keeping the initial flights.
        best = schedules.timetable(routes, [p["assign"] for p in initial], fixed=starts)
        objective = validate(instance, best)
    elif not any(r["parcels"] for r in routes):
        best = baselines
        objective = validate(instance, best)

    if objective == common_bound:
        return result("optimal", common_bound, "common")
    if progress is not None:
        progress(objective, common_bound)
    if perf_counter() >= deadline:
        return result("limit", common_bound)

    caps = [horizon - r["base"] for r in routes]
    if objective is not None:
        caps = [min(cap, int(objective) - common_bound + lb)
                for cap, lb in zip(caps, lower)]

    outgoing, incoming = defaultdict(list), defaultdict(list)
    build_start = perf_counter()
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", 0)
        env.start()
        with gp.Model("liu_fixed_fstsp", env=env) as model:
            model.ModelSense = GRB.MINIMIZE
            model.Params.Threads = threads
            model.Params.Seed = 13
            model.Params.MIPGap = 0
            model.Params.MIPGapAbs = 0
            model.Params.FeasibilityTol = 1e-8
            model.Params.IntFeasTol = 1e-8
            model.Params.OptimalityTol = 1e-8
            teams = []
            for k, (route, b, tail, last, lb, ub) in enumerate(
                    zip(routes, baselines, tails, latest, lower, caps)):
                team = _add_team(
                    model, route, b, tail, last, lb, ub, instance, outgoing, incoming,
                    best[k] if best is not None else None, deadline, stats,
                )
                if team is None:
                    stats["build_seconds"] = perf_counter() - build_start
                    return result("limit", common_bound)
                teams.append(team)

            complete = _add_fleet_flow(
                model, outgoing, incoming, horizon, instance["M"], deadline, stats
            )
            if complete:
                model.update()
            stats["build_seconds"] = perf_counter() - build_start
            if not complete or perf_counter() >= deadline:
                return result("limit", common_bound)

            stats.update(variables=model.NumVars, binaries=model.NumBinVars,
                         constraints=model.NumConstrs)

            def callback(m, where):
                if where == GRB.Callback.MIP:
                    incumbent = m.cbGet(GRB.Callback.MIP_OBJBST)
                    bound = m.cbGet(GRB.Callback.MIP_OBJBND)
                    progress(incumbent if _finite_bound(incumbent) else objective,
                             max(common_bound, bound) if _finite_bound(bound)
                             else common_bound)

            optimize_start = perf_counter()
            model.Params.TimeLimit = max(0, deadline - optimize_start)
            if progress is None:
                model.optimize()
            else:
                model.optimize(callback)
            stats.update(optimizer_seconds=perf_counter() - optimize_start,
                         solver_status=model.Status, nodes=model.NodeCount)
            if model.Status in (GRB.INFEASIBLE, GRB.INF_OR_UNBD):
                if best is not None:
                    raise ArithmeticError("Network infeasibility contradicts the initial schedule")
                return result("infeasible", None)
            if model.SolCount:
                candidate = _witness(instance, model, teams)
                value = validate(instance, candidate)
                if value > model.ObjVal + TOL:
                    raise ArithmeticError("Physical replay exceeds the network objective")
                if objective is None or value < objective:
                    objective, best = value, candidate

            bound = common_bound
            if _finite_bound(model.ObjBound):
                bound = max(bound, model.ObjBound)
            if objective is not None:
                if bound > objective + TOL:
                    raise ArithmeticError("Network lower bound exceeds the verified schedule")
                if ceil(bound - TOL) >= objective:
                    return result("optimal", objective)
            return result("limit", bound)
