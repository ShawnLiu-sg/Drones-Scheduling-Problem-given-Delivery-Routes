"""Lower bounds from mandatory replenishment energy at alternative route delays.

For a route allowed to finish at base+d, every parcel induces a collection of
possible (stop, earliest start, latest start, duration) windows. If the start
can range from f to l, the smallest possible overlap with [a,b) is
max(0, min(b-a, duration, b-l, f+duration-a)). Taking the minimum over stops
is a per-parcel compulsory energy. Distinct forced stops use distinct flights.
Parcel-weighted energy also bounds flight energy because a flight carries at
most Q. A multiple-choice master selects one delay option per route subject
to the fleet's M*(b-a) supply in every interval. Any physical schedule induces
a feasible master choice of no greater cost, including the overflow option.
"""
from collections import defaultdict
from math import ceil, isfinite
from pathlib import Path
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB

from scheduling.study.reference import _latest, validate


def _closure(values, arcs):
    for _ in values:
        changed = False
        for i, j, lag in arcs:
            if values[i] + lag > values[j]:
                values[j] = values[i] + lag
                changed = True
        if not changed:
            return values
    raise ValueError("Inconsistent fixed-route time network")


def _event_envelope(route):
    n, s = len(route["tau"]), len(route["sorties"])
    a = list(range(n))
    d = list(range(n, 2*n))
    l = list(range(2*n, 2*n+s))
    u = list(range(2*n+s, 2*n+2*s))
    arcs = []
    for j in range(1,n):
        travel = route["travel"][j-1]
        arcs.extend(((d[j-1],a[j],travel), (a[j],d[j-1],-travel),
                     (a[j],d[j],route["service"][j])))
    for j, f in enumerate(route["sorties"]):
        for event, stop in ((l[j],f["l"]), (u[j],f["b"])):
            arcs.extend(((a[stop],event,0), (event,d[stop],0)))
        arcs.extend(((l[j],u[j],f["F"]), (u[j],l[j],-f["E"])))
        if j:
            arcs.append((u[j-1],l[j],0))
    slots = 2*n + 2*s
    early = [0]*slots
    radially_dominated = all(length >= abs(x-y) for length,x,y in
                             zip(route["travel"],route["tau"],route["tau"][1:]))
    if radially_dominated:
        for p in route["parcels"]:
            dest = _latest(route,p)
            event = l[p["dest"]] if p["kind"]=="sortie" else d[dest]
            service = route["service"][dest] if p["kind"]=="truck" else 0
            early[event] = max(early[event],p["r"]+route["tau"][dest]+service)
    early = _closure(early,arcs)
    tail = [float("-inf")]*slots
    tail[d[-1]] = 0
    tail = _closure(tail,[(j,i,lag) for i,j,lag in arcs])
    return early,tail,a,d,l


def _windows(route, horizon, capacity):
    early,tail,a,d,l = _event_envelope(route)
    bounds = []
    for parcel in route["parcels"]:
        eligible = {}
        for i in range(1, min(_latest(route,parcel),len(a)-2)+1):
            tau = route["tau"][i]
            if not tau or parcel["w"] > capacity:
                continue
            first = max(0,parcel["r"],early[a[i]]-tau)
            after = tail[d[i]] + tau
            if parcel["kind"]=="sortie":
                after = max(after,tail[l[parcel["dest"]]]+tau)
            elif parcel["dest"]==i:
                after += route["service"][i]
            eligible[i] = (first, route["base"]-after, horizon-2*tau)
        bounds.append(eligible)
    return bounds, max(0,early[d[-1]]-route["base"])


def _classify(route, bounds, delay, capacity):
    weights = [p["w"] for p in route["parcels"]]
    sites = [{i: (first, min(end+delay,latest)) for i,(first,end,latest) in row.items()
              if first <= min(end+delay,latest)} for row in bounds]
    while True:
        if any(not opt for opt in sites):
            return None
        mandatory = defaultdict(list)
        for j, opt in enumerate(sites):
            if len(opt)==1:
                mandatory[next(iter(opt))].append(j)
        for i, jobs in mandatory.items():
            load = sum(weights[j] for j in jobs)
            first = max(sites[j][i][0] for j in jobs)
            last = min(sites[j][i][1] for j in jobs)
            if load > capacity or first > last:
                return None
            for j, opt in enumerate(sites):
                if i in opt:
                    left, right = opt[i]
                    left, right = max(left,first), min(right,last)
                    if left > right or (j not in jobs and weights[j]+load > capacity):
                        del opt[i]
                    else:
                        opt[i] = (left,right)
        if sum(map(len,mandatory.values())) == sum(len(opt)==1 for opt in sites):
            if any(not opt for opt in sites):
                return None
            break

    # Each tuple gives the envelope of one parcel's possible flights.  A
    # forced stop may hold several parcels but contributes one distinct trip.
    forced = [(max(sites[j][i][1] for j in jobs),
               min(sites[j][i][0] for j in jobs)+2*route["tau"][i],
               2*route["tau"][i]) for i,jobs in mandatory.items()]
    pieces = []
    for j,opt in enumerate(sites):
        latest = max(r for _,r in opt.values())
        finish = min(f+2*route["tau"][i] for i,(f,_) in opt.items())
        duration = min(2*route["tau"][i] for i in opt)
        pieces.append((latest,finish,duration,weights[j],mandatory.keys().isdisjoint(opt)))
    return forced,pieces


def _energy(profile, intervals, capacity):
    forced,pieces = profile
    energies = []
    for a,b in intervals:
        width=b-a
        fixed=sum(max(0,min(width,p,b-l,e-a)) for l,e,p in forced)
        volume=free_volume=free_peak=peak=0
        for l,e,p,w,free in pieces:
            part=max(0,min(width,p,b-l,e-a))
            volume += w*part
            peak=max(peak,part)
            if free:
                free_volume+=w*part
                free_peak=max(free_peak,part)
        energies.append(max(ceil(volume/capacity-1e-9),peak,
                            fixed+max(free_peak,ceil(free_volume/capacity-1e-9))))
    return tuple(energies)



def _master_bound(tables, capacities, deadline, env, threads):
    """Solve a relaxation on the supplied teams with the full fleet capacity."""
    baseline=sum(table[0][0] for table in tables)
    if perf_counter()>=deadline:
        return baseline
    active=[i for i,rhs in enumerate(capacities)
            if sum(table[0][1][i] for table in tables)>rhs]
    if not active:
        return baseline
    with gp.Model("mandatory_flight_energy",env=env) as m:
        m.Params.Threads=threads
        m.Params.MIPGap=0
        m.Params.MIPGapAbs=0
        rows=[[] for _ in active]
        for table in tables:
            choices=[]
            for cost,amounts in table:
                var=m.addVar(vtype=GRB.BINARY,obj=cost)
                choices.append(var)
                for row,i in zip(rows,active):
                    if amounts[i]:
                        row.append((amounts[i],var))
            m.addLConstr(gp.quicksum(choices)==1)
        for row,i in zip(rows,active):
            m.addLConstr(gp.quicksum(w*x for w,x in row)<=capacities[i])
        m.Params.TimeLimit=max(0,deadline-perf_counter())
        m.optimize()
        if m.Status!=GRB.LOADED and isfinite(m.ObjBound) and abs(m.ObjBound)<1e90:
            return max(baseline,ceil(m.ObjBound-1e-6))
        return baseline


def energy_cuts(instance,time_limit=2.0,max_extra=12,threads=1,*,
                joint_bound=False,other_bounds=False):
    """Return valid interval cuts and delay costs for integration in a master.

    Select one (min_delay, energy_vector) per team and require min_delay <=
    actual route delay. Each vector coordinate sums to at most M*(b-a).
    A final zero-energy choice covers delays above the tabulated range.
    Bound-only timeouts retain the route lower bounds already proved.
    """
    start=perf_counter()
    deadline=start+time_limit
    routes=instance["routes"]
    horizon,fleet,capacity=(instance[key] for key in ("H","M","Q"))
    windows=[_windows(r,horizon,capacity) for r in routes]
    local=[bound for _,bound in windows]
    last=min(horizon,max((r["base"]+v+max_extra+max(r["tau"])
                          for r,v in zip(routes,local)),default=0))
    intervals=[(a,a+w) for w in (1,4,8,16) for a in range(0,last-w+1,max(1,w//4))]
    neutral=(0,)*len(intervals)
    tables=[[(lb,neutral)] for lb in local]
    for k,(route,(raw,low)) in enumerate(zip(routes,windows)):
        if perf_counter()>=deadline:
            break
        options=[]
        previous=None
        for delay in range(low,low+max_extra+1):
            if perf_counter()>=deadline:
                break
            profile=_classify(route,raw,delay,capacity)
            if profile is not None and profile!=previous:
                vector=_energy(profile,intervals,capacity)
                if not options or vector!=options[-1][1]:
                    options.append((delay,vector))
            previous=profile
        else:
            options.append((low+max_extra+1,neutral))
            tables[k]=options
            local[k]=options[0][0]
    keep=[i for i,(a,b) in enumerate(intervals)
          if sum(table[0][1][i] for table in tables)>fleet*(b-a)]
    intervals=[intervals[i] for i in keep]
    capacities=[fleet*(b-a) for a,b in intervals]
    filtered=[]
    for table in tables:
        levels=[]
        for cost,vec in table:
            values=tuple(vec[i] for i in keep)
            if not levels or levels[-1][1] != values:
                levels.append((cost,values))
        filtered.append(levels)
    bound=sum(local)
    others=[bound-lb for lb in local]
    if keep and (joint_bound or other_bounds) and perf_counter()<deadline:
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag",0)
            env.start()
            if joint_bound:
                bound=_master_bound(filtered,capacities,deadline,env,threads)
            if other_bounds:
                for k in range(len(routes)):
                    if perf_counter()>=deadline:
                        break
                    others[k]=_master_bound(filtered[:k]+filtered[k+1:],capacities,
                                            deadline,env,threads)
    return dict(intervals=intervals, options=filtered, route_bounds=local,
                other_bounds=others, bound=bound,
                stats=dict(energy_rows=len(keep),choices=sum(map(len,filtered)),
                           seconds=perf_counter()-start))


def mandatory_bound(instance,time_limit=2.0,max_extra=12,threads=1):
    """Return a certified integer LB, route-wise LBs, and the common witness."""
    start=perf_counter()
    warm=instance.get("warm")
    ub=validate(instance,warm) if warm is not None else None
    cert=energy_cuts(instance,max(0,time_limit-(perf_counter()-start)),max_extra,
                     threads,joint_bound=True)
    bound=cert["bound"]
    return dict(status="optimal" if ub==bound else "limit",objective=ub,bound=bound,
                routes=warm,route_bounds=cert["route_bounds"],origin="algorithm",
                stats=cert["stats"]|{"seconds":perf_counter()-start})
