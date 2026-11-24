"""Deterministic, paired synthetic fixed-route benchmarks, including a common incumbent."""
import argparse
from itertools import product
import json
from pathlib import Path
from random import Random
from .model import baseline, eligible, timetable

def route(rng, stops, k, seed):
    radius = [0] + [rng.randint(2, 5) for _ in range(stops)] + [0]
    service = [0] + [rng.randint(1, 2) for _ in range(stops)] + [0]
    travel = [abs(a-b)+rng.randint(1, 3) for a,b in zip(radius, radius[1:])]
    legs = [(1, 3)] if stops < 5 else [(1, 3), (3, 5)]
    sorties = []
    for l,b in legs:
        road = sum(travel[l:b]) + sum(service[l+1:b])
        flying = max(abs(radius[l]-radius[b])+1, road+rng.randint(-2, 1))
        sorties.append(dict(l=l,b=b,F=flying,E=max(flying,road)+rng.randint(0, 2)))
    parcels = [dict(kind="sortie",dest=0,r=rng.randint(2, 5),w=rng.randint(1,2)),
               dict(kind="truck",dest=stops,r=rng.randint(7, 12),w=rng.randint(1,2))]
    if stops>=4:
        parcels.append(dict(kind="truck",dest=2,r=rng.randint(4, 8),w=rng.randint(1,2)))
    if stops>=5:
        parcels.append(dict(kind="sortie",dest=1,r=rng.randint(6,11),w=rng.randint(1,2)))
    rt=dict(travel=travel,service=service,tau=radius,sorties=sorties,parcels=parcels)
    rt["base"]=baseline(rt)
    return rt


def warm_schedule(routes, m, capacity, seed):
    rng=Random(81041+seed)
    pools=[]
    for rt in routes:
        poss=[v for v in product(*(eligible(rt,p) for p in rt["parcels"]))
              if all(sum(p["w"] for p,s in zip(rt["parcels"],v) if s==i)<=capacity
                     for i in set(v))]
        pools.append(sorted(poss,key=lambda v:(len(set(v)),sum(rt["tau"][i] for i in set(v)),v)))
    best, best_c=None, float("inf")
    for it in range(80):
        assignments=[ps[0 if it==0 else rng.randrange(len(ps))] for ps in pools]
        # Each route's meetings stay in visit order; interleave routes and assign free lanes.
        jobs=[sorted(set(v)) for v in assignments]
        lanes=[[] for _ in range(m)]
        clocks=[0]*m
        while any(jobs):
            active=[k for k,q in enumerate(jobs) if q]
            k=rng.choice(active)
            i=jobs[k].pop(0)
            rt=routes[k]
            lane=min(range(m),key=lambda d:clocks[d])
            ready=max(p["r"] for p,s in zip(rt["parcels"],assignments[k]) if s==i)
            lanes[lane].append((k,i))
            clocks[lane]=max(clocks[lane],ready)+2*rt["tau"][i]
        sol=timetable(routes,assignments,chains=lanes)
        if sol is not None:
            cost=sum(w["D"][-1]-rt["base"] for rt,w in zip(routes,sol))
            if cost<best_c:
                best,best_c=sol,cost
    if best is None:
        raise ValueError("No common feasible starter found")
    return best


def instance(k, stops, m, seed):
    rng=Random(20261118 + 10000*k+100*stops+seed)
    routes=[route(rng,stops,j,seed) for j in range(k)]
    warm=warm_schedule(routes,m,3,seed)
    ub=sum(s["D"][-1]-r["base"] for r,s in zip(routes,warm))
    return dict(name=f"k{k}-p{stops}-m{m}-s{seed}",M=m,Q=3,
                H=max(r["base"] for r in routes)+ub,routes=routes,warm=warm,
                parameters=dict(teams=k,stops=stops,seed=seed,initial_delay=ub))
