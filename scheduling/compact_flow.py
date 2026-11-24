"""Exact cumulative-threshold scheduling for fixed FSTSP routes and a shared fleet."""

from collections import defaultdict
from math import ceil
from pathlib import Path
import sys
from time import perf_counter

from scheduling import liu_flow as liu
from gurobipy import GRB
import gurobipy as gp


class Clock:
    """Survival function of an integer event on [lo,hi]."""
    def __init__(self, model, lo, hi):
        self.lo, self.hi = lo, hi
        self.tail = list(model.addVars(hi - lo, ub=1).values())
        for left, right in zip(self.tail, self.tail[1:]):
            model.addLConstr(right <= left)

    def at(self, time):
        if time <= self.lo:
            return 1
        if time > self.hi:
            return 0
        return self.tail[time - self.lo - 1]


class Dispatch:
    """Optional integral departure encoded by descending binary indicators."""
    def __init__(self, model, first, last, tau):
        self.first, self.last, self.tau = first, last, tau
        self.tail = list(model.addVars(last - first + 1, vtype=GRB.BINARY).values())
        self.used = self.tail[0]
        for left, right in zip(self.tail, self.tail[1:]):
            model.addLConstr(right <= left)

    def at(self, tick):
        if tick > self.last:
            return 0
        return self.tail[max(0, tick - self.first)]


def _clock_graph(route, baseline, cap, horizon):
    end = len(route["service"]) - 1
    lo = {("D", i): value for i, value in enumerate(baseline["D"])}
    arcs = [(("D", i-1), ("D", i), route["travel"][i-1] + route["service"][i])
            for i in range(1, end+1)]
    for j, sortie in enumerate(route["sorties"]):
        for tag, stop in (("L", sortie["l"]), ("U", sortie["b"])):
            event = (tag, j)
            lo[event] = baseline[tag][j]
            if stop:
                arcs.append((("D", stop-1), event, route["travel"][stop-1]))
            arcs.append((event, ("D", stop), 0))
        arcs.extend([(("L", j), ("U", j), sortie["F"]),
                     (("U", j), ("L", j), -sortie["E"])])
        if j:
            arcs.append((("U", j-1), ("L", j), 0))
    hi = {e:horizon for e in lo}
    hi["D", 0] = 0
    hi["D", end] = min(horizon, route["base"] + cap)
    for _ in lo:
        changed = False
        for a, b, lag in reversed(arcs):
            if hi[a] > hi[b]-lag:
                hi[a] = hi[b]-lag
                changed = True
        if not changed:
            break
    return lo, hi, arcs


def _team(model, inst, rt, base, latest, cap, lower, deadline):
    lo, hi, arcs = _clock_graph(rt, base, cap, inst["H"])
    if any(lo[e] > hi[e] for e in lo):
        return None
    clocks = {e:Clock(model, lo[e], hi[e]) for e in lo}
    for a, b, lag in arcs:
        left, right = clocks[a], clocks[b]
        for t in range(max(left.lo+1, right.lo-lag+1), left.hi+1):
            model.addLConstr(left.at(t) <= right.at(t+lag))

    choices = [[] for _ in latest]
    stops = {}
    for i in range(1, max(latest, default=0)+1):
        if perf_counter() >= deadline:
            raise TimeoutError
        tau = rt["tau"][i]
        if tau <= 0:
            continue
        before, depart = clocks["D", i-1], clocks["D", i]
        road = rt["travel"][i-1]
        last = min(depart.hi-tau, inst["H"]-2*tau)
        ps = rt["parcels"]
        eligible = [p for p, limit in enumerate(latest) if i <= limit and ps[p]["r"] <= last]
        if not eligible:
            continue
        first = max(0, before.lo+road-tau, min(ps[p]["r"] for p in eligible))
        if first > last:
            continue
        h = Dispatch(model, first, last, tau)
        stops[i] = h
        # The truck must have arrived by receipt, and can only leave afterwards.
        for t in range(max(before.lo+1, first-road+tau+1), before.hi+1):
            model.addLConstr(before.at(t) <= 1-h.used + h.at(t+road-tau))
        for t in range(max(first, depart.lo-tau+1), last+1):
            model.addLConstr(depart.at(t+tau) >= h.at(t))
        load = model.addVars(eligible, vtype=GRB.BINARY)
        model.addLConstr(load.sum() >= h.used)
        for p in eligible:
            parcel = ps[p]
            choices[p].append((i,load[p]))
            model.addLConstr(load[p] <= h.at(parcel["r"]))
            event, lag = None, tau
            if parcel["kind"] == "truck" and parcel["dest"] == i and rt["service"][i]:
                event, lag = depart, tau + rt["service"][i]
            if parcel["kind"] == "sortie" and rt["sorties"][parcel["dest"]]["l"] == i:
                event = clocks["L", parcel["dest"]]
            # When the parcel is delivered further along the truck route, the
            # truck precedence constraints already enforce its receipt by launch.
            if event is not None:
                for t in range(max(first,event.lo-lag+1), last+1):
                    model.addLConstr(h.at(t) + load[p] - 1 <= event.at(t+lag))
        for release in sorted({max(first,ps[p]["r"]) for p in eligible}):
            subset = [p for p in eligible if max(first,ps[p]["r"]) >= release]
            model.addLConstr(gp.quicksum(ps[p]["w"]*load[p] for p in subset)
                             <= inst["Q"]*h.at(release))

    for options in choices:
        model.addLConstr(gp.quicksum(q for _,q in options) == 1)
    end = clocks["D",len(rt["service"])-1]
    delay = end.lo - rt["base"] + gp.quicksum(end.tail)
    model.addLConstr(delay >= lower)
    return stops, choices, delay, end


def _trial(env, inst, setup, budget, deadline, threads, energy=None):
    baselines, latest, lower = setup
    with gp.Model("cumulative_thresholds", env=env) as m:
        m.Params.Threads = threads
        m.Params.Seed = 13
        m.Params.MIPGap = 0
        m.Params.MIPGapAbs = 0
        m.Params.IntFeasTol = 1e-8
        m.Params.FeasibilityTol = 1e-8
        teams = []
        for k, rt in enumerate(inst["routes"]):
            cap = budget-sum(lower)+lower[k]
            if energy is not None:
                cap = min(cap, budget-energy["other_bounds"][k])
            row = _team(m,inst,rt,baselines[k],latest[k],cap,lower[k],deadline)
            if row is None:
                return "empty", None, budget+1, {}
            teams.append(row)
        fleet_rows = defaultdict(gp.LinExpr)
        for stops, _, _, _ in teams:
            for h in stops.values():
                for t in range(h.first,h.last+2*h.tau):
                    fleet_rows[t] += h.at(t-2*h.tau+1)-h.at(t+1)
        for load in fleet_rows.values():
            m.addLConstr(load <= inst["M"])
        if energy is not None:
            for j, (a, b) in enumerate(energy["intervals"]):
                occupancy = gp.LinExpr()
                for rt, team, options in zip(inst["routes"], teams, energy["options"]):
                    occupancy += options[0][1][j]
                    for (_, old), (delay_level, new) in zip(options, options[1:]):
                        delta = new[j] - old[j]
                        if delta:
                            occupancy += delta * team[3].at(rt["base"] + delay_level)
                m.addLConstr(occupancy <= inst["M"] * (b - a))
        delay = gp.quicksum(row[2] for row in teams)
        m.addLConstr(delay <= budget)
        m.setObjective(delay)
        remaining = deadline - perf_counter()
        if remaining <= 0:
            raise TimeoutError
        m.Params.TimeLimit = remaining
        m.optimize()
        statistics = dict(variables=m.NumVars, binary=m.NumBinVars, rows=m.NumConstrs,
                          nodes=m.NodeCount, runtime=m.Runtime)
        if m.Status == GRB.INFEASIBLE:
            return "empty",None,budget+1,statistics
        bound = min(budget+1, ceil(m.ObjBound-1e-6)) if liu._finite_bound(m.ObjBound) else sum(lower)
        witness = None
        if m.SolCount:
            fixed, batches = [], []
            for stops, loads, _, _ in teams:
                starts = {i:h.first+sum(v.X>0.5 for v in h.tail)-1
                          for i,h in stops.items() if h.used.X>0.5}
                fixed.append(starts)
                batches.append([next(i for i,v in options if v.X>0.5) for options in loads])
            witness = liu.schedules.timetable(inst["routes"],batches,fixed=fixed)
            ub = liu.validate(inst,witness)
            if ub > m.ObjVal+1e-5:
                raise ArithmeticError(("Inconsistent threshold witness",ub,m.ObjVal))
        return "solved" if m.Status == GRB.OPTIMAL else "limit",witness,bound,statistics


def solve(instance,time_limit=120,threads=1,*,progress=None):
    """Return a verified schedule and certified integer delay lower bound."""
    started = perf_counter()
    deadline = started + time_limit
    routes = instance["routes"]
    if "H" not in instance:
        return liu.solve(instance,time_limit,threads,progress=progress)
    best = instance.get("warm")
    upper = liu.validate(instance,best) if best is not None else None
    lower = 0
    trials = []

    def report(status):
        if progress:
            progress(upper, lower)
        return dict(status=status,objective=upper,bound=None if status=="infeasible" else lower,
                    routes=best,stats=dict(trials=trials,seconds=perf_counter()-started),
                    origin="algorithm")
    bases = liu.schedules.timetable(routes,[[] for _ in routes])
    if bases is None:
        return report("infeasible")
    latest = [[min(liu._latest(r,p),len(r["tau"])-2) for p in r["parcels"]] for r in routes]
    lbs = [liu._team_lower_bound(r,b,liu._road_tail(r),p,instance["Q"])
           for r,b,p in zip(routes,bases,latest)]
    if any(v is None or r["base"]+v>instance["H"] for v,r in zip(lbs,routes)):
        return report("infeasible")
    lower = sum(lbs)
    if best is None and not any(r["parcels"] for r in routes):
        best, upper = bases, liu.validate(instance,bases)
    if upper == lower:
        return report("optimal")
    if progress:
        progress(upper,lower)
    highest = upper-1 if upper is not None else sum(instance["H"]-r["base"] for r in routes)
    budget = int(min(highest, lower+max(1,(highest-lower+1)//8)))
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag",0)
        env.start()
        while perf_counter() < deadline:
            try:
                status,witness,lb,stats = _trial(env,instance,(bases,latest,lbs),budget,deadline,threads)
            except TimeoutError:
                break
            trials.append(dict(budget=budget,status=status,**stats))
            if witness is not None:
                value = liu.validate(instance,witness)
                if upper is None or upper>value:
                    best,upper = witness,value
            lower = max(lower,lb)
            if upper is not None and lower >= upper:
                if lower != upper:
                    raise ArithmeticError(("Certificate conflict",lower,upper))
                return report("optimal")
            if status == "empty" and budget == highest:
                return report("infeasible")
            if progress:
                progress(upper,lower)
            if status != "empty":
                break
            budget = min(highest, max(lower, 2*budget-sum(lbs)+2))
    return report("limit")
