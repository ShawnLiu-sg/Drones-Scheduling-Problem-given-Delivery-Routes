"""Pilot MILP using one potential replenishment flight per late parcel.

A batch is identified by its lowest-index parcel. The selected flights form
at most M directed paths. Positive round-trip durations exclude cycles, so
the paths can be assigned to M distinct replenishment drones.
"""
from collections import defaultdict
from pathlib import Path
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB

from scheduling.study import model as clock
from scheduling.study import reference


def _route_events(model, route, endcap):
    n = len(route["service"])
    order = range(n)
    clean_route = dict(route, parcels=[])
    baseline = clock.timetable([clean_route], [[]])[0]
    events = {s: model.addVars(len(baseline[s]), lb=baseline[s], ub=endcap, name=s)
              for s in ("A", "D", "L", "U")}
    a, d, l, u = [events[s] for s in ("A", "D", "L", "U")]
    model.addConstr(a[0] == 0)
    model.addConstr(d[0] == 0)
    for i in order[1:]:
        model.addConstr(a[i] == d[i-1] + route["travel"][i-1])
        model.addConstr(d[i] >= a[i] + route["service"][i])
    for j, sortie in enumerate(route["sorties"]):
        launch, recover = sortie["l"], sortie["b"]
        model.addConstr(l[j] >= a[launch])
        model.addConstr(l[j] <= d[launch])
        model.addConstr(u[j] >= a[recover])
        model.addConstr(u[j] <= d[recover])
        model.addConstr(u[j] - l[j] >= sortie["F"])
        model.addConstr(u[j] - l[j] <= sortie["E"])
        if j:
            model.addConstr(l[j] >= u[j-1])
    return events, baseline


def _create_team(model, route, endcap, H, capacity, fixed):
    events, baseline = _route_events(model, route, endcap)
    a, d, l = [events[s] for s in ("A", "D", "L")]
    parcels = route["parcels"]
    count = len(parcels)
    n = len(a)
    last = [reference._latest(route, p) for p in parcels]
    tail = [0]*n
    for i in range(n-2, -1, -1):
        tail[i] = tail[i+1] + route["travel"][i] + route["service"][i+1]
    slots = []
    for j, p in enumerate(parcels):
        possible = [i for i in range(1, last[j]+1)
                    if route["tau"][i] > 0
                    and p["r"] + 2*route["tau"][i] <= H
                    and p["r"] + route["tau"][i] + tail[i] <= endcap
                    and (fixed is None or fixed[j] == i)]
        slots.append(set(possible))
    active = model.addVars(count, vtype=GRB.BINARY, name="batch")
    station = [model.addVars(sorted(s), vtype=GRB.BINARY, name="station") for s in slots]
    members = model.addVars([(p,b) for p in range(count) for b in range(p+1)
                             if p == b or slots[p] & slots[b]], vtype=GRB.BINARY, name="member")
    start = model.addVars(count, ub=H, name="start")
    flight = model.addVars(count, ub=H, name="one_way")
    meet = model.addVars(count, ub=endcap, name="handoff")
    for i in sorted(set().union(*slots)):
        model.addConstr(gp.quicksum(station[b][i] for b in range(count) if i in slots[b]) <= 1)
    for p, parcel in enumerate(parcels):
        model.addConstr(members.sum(p, "*") == 1)
        for b in range(p+1):
            if (p,b) not in members:
                continue
            z = members[p,b]
            if p == b:
                model.addConstr(z == active[b])
            else:
                model.addConstr(z <= gp.quicksum(station[b][i] for i in sorted(slots[p] & slots[b])))
            model.addConstr(start[b] >= parcel["r"] * z)
            if parcel["kind"] == "sortie":
                model.addGenConstrIndicator(z, True, l[parcel["dest"]] >= meet[b])
            else:
                dest = parcel["dest"]
                # If handed over earlier on the route, this follows from truck travel and service.
                model.addGenConstrIndicator(z, True, d[dest] >= meet[b] + route["service"][dest])
    for b in range(count):
        choose = station[b]
        model.addConstr(choose.sum() == active[b])
        model.addConstr(gp.quicksum(parcels[p]["w"]*members[p,b] for p in range(b,count)
                                    if (p,b) in members) <= capacity*active[b])
        model.addConstr(flight[b] == gp.quicksum(route["tau"][i]*v for i,v in choose.items()))
        model.addConstr(meet[b] == start[b] + flight[b])
        model.addConstr(meet[b] + flight[b] <= H*active[b])
        if slots[b]:
            model.addConstr(meet[b] >= gp.quicksum(baseline["A"][i]*v for i,v in choose.items()))
            model.addConstr(d[n-1] >= meet[b] + gp.quicksum(tail[i]*v for i,v in choose.items()))
        for i,v in choose.items():
            model.addGenConstrIndicator(v, True, a[i] <= meet[b])
            model.addGenConstrIndicator(v, True, d[i] >= meet[b])
    return dict(events=events, active=active, station=station, members=members, start=start,
                flight=flight, meet=meet, terminal=d[n-1])


def _warm_start(instance, teams, roots, edges):
    intervals = []
    for k,(team, rt, w) in enumerate(zip(teams,instance["routes"],instance["warm"])):
        for key,variables in team["events"].items():
            for j,v in variables.items():
                v.Start = w[key][j]
        assignment = w["assign"]
        starts = {int(i):v for i,v in w["starts"].items()}
        leaders = {i: assignment.index(i) for i in set(assignment)}
        for b,i in enumerate(assignment):
            used = b == leaders[i]
            t = starts[i] if used else 0
            tau = rt["tau"][i] if used else 0
            team["active"][b].Start = used
            team["start"][b].Start = t
            team["flight"][b].Start = tau
            team["meet"][b].Start = t + tau
            for h,v in team["station"][b].items():
                v.Start = int(used and h == i)
            if used:
                intervals.append((t,t+2*tau,(k,b)))
        for (p,b),v in team["members"].items():
            v.Start = int(leaders[assignment[p]] == b)
    for v in list(roots.values()) + list(edges.values()):
        v.Start = 0
    lanes = []
    for t,end,job in sorted(intervals):
        available = next((j for j,(free,_) in enumerate(lanes) if free <= t), None)
        if available is None:
            roots[job].Start = 1
            lanes.append((end,job))
        else:
            edges[lanes[available][1],job].Start = 1
            lanes[available] = (end,job)


def run(instance, seconds, fixed_assignments=None):
    begin = perf_counter()
    warm = instance.get("warm")
    if warm is not None and fixed_assignments is not None:
        if [list(r["assign"]) for r in warm] != [list(row) for row in fixed_assignments]:
            warm = None
    incumbent = reference.validate(instance,warm) if warm is not None else None
    H = instance["H"]
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag",0)
        env.start()
        with gp.Model("pilot_batch_path_cover",env=env) as mip:
            mip.Params.Threads = 1
            mip.Params.Seed = 13
            mip.Params.MIPGap = mip.Params.MIPGapAbs = 0
            mip.Params.IntFeasTol = mip.Params.FeasibilityTol = 1e-8
            teams = []
            for k,r in enumerate(instance["routes"]):
                endcap = H if incumbent is None else min(H, r["base"] + incumbent)
                teams.append(_create_team(mip, r, endcap, H, instance["Q"],
                             fixed_assignments[k] if fixed_assignments is not None else None))
            jobs = [(k,b) for k,team in enumerate(teams) for b in team["active"]]
            roots = {a: mip.addVar(vtype=GRB.BINARY, name="root") for a in jobs}
            edges = {}
            out_arcs, in_arcs = defaultdict(list), defaultdict(list)
            for a in jobs:
                k,b = a
                for c in jobs:
                    if a == c:
                        continue
                    j,s = c
                    v = mip.addVar(vtype=GRB.BINARY,name="successor")
                    edges[a,c] = v
                    out_arcs[a].append(v)
                    in_arcs[c].append(v)
                    mip.addGenConstrIndicator(v,True,teams[j]["start"][s] >=
                                              teams[k]["meet"][b]+teams[k]["flight"][b])
            for a in jobs:
                active = teams[a[0]]["active"][a[1]]
                mip.addConstr(gp.quicksum(out_arcs[a]) <= active)
                mip.addConstr(roots[a] + gp.quicksum(in_arcs[a]) == active)
            mip.addConstr(gp.quicksum(roots.values()) <= instance["M"])
            mip.setObjective(gp.quicksum(team["terminal"]-r["base"]
                                        for team,r in zip(teams,instance["routes"])))
            if warm is not None:
                _warm_start(instance, teams, roots, edges)
            mip.Params.TimeLimit = max(0, seconds - (perf_counter()-begin))
            mip.optimize()
            status = {GRB.OPTIMAL:"optimal", GRB.INFEASIBLE:"infeasible"}.get(mip.Status,"limit")
            plan, obj = warm, incumbent
            if mip.SolCount:
                assignments = []
                sites = {}
                for k,team in enumerate(teams):
                    active_sites = {b:next(i for i,v in choices.items() if v.X > .5)
                                    for b,choices in enumerate(team["station"]) if team["active"][b].X > .5}
                    sites.update({(k,b):(k,i) for b,i in active_sites.items()})
                    row = [0]*len(team["station"])
                    for (p,b),v in team["members"].items():
                        if v.X > .5:
                            row[p] = active_sites[b]
                    assignments.append(row)
                successor = {a:c for (a,c),v in edges.items() if v.X > .5}
                lanes = []
                for head,v in roots.items():
                    if v.X > .5:
                        lane = [sites[head]]
                        while head in successor:
                            head = successor[head]
                            lane.append(sites[head])
                        lanes.append(lane)
                schedule = clock.timetable(instance["routes"],assignments,chains=lanes)
                val = reference.validate(instance,schedule)
                if val > mip.ObjVal + 1e-5 or (status == "optimal" and abs(val-mip.ObjVal)>1e-5):
                    raise ArithmeticError(("solver objective",mip.ObjVal,"physical schedule",val))
                if obj is None or val < obj:
                    plan,obj = schedule,val
            if status == "infeasible" and warm is not None:
                raise ArithmeticError("Infeasible certificate conflicts with replayed warm schedule")
            bound = None if status == "infeasible" else max(0.,mip.ObjBound)
            if status=="optimal":
                bound = obj
            if obj is not None and bound is not None and bound > obj+1e-5:
                raise ArithmeticError(("inconsistent bounds",bound,obj))
            return dict(status=status,objective=obj,bound=bound,routes=plan,
                        stats=dict(seconds=perf_counter()-begin,nodes=mip.NodeCount,variables=mip.NumVars,
                                   binaries=mip.NumBinVars,constraints=mip.NumConstrs,
                                   indicators=mip.NumGenConstrs,solver_status=mip.Status))
