"""Exhaustive integer-grid check on the two-team cases, with no optimizer calls."""
from itertools import product
import json
from pathlib import Path
from .model import eligible, realize

def local_profiles(r, Q, horizon, M):
    profiles={}
    for assign in product(*(eligible(r,p) for p in r["parcels"])):
        sites=sorted(set(assign))
        bundles=[[p for p,i in zip(r["parcels"],assign) if i==s] for s in sites]
        if any(sum(p["w"] for p in ps)>Q for ps in bundles):
            continue
        domains=[range(max(p["r"] for p in ps),horizon-2*r["tau"][s]+1)
                 for s,ps in zip(sites,bundles)]
        for times in product(*domains):
            masks=[0]*M
            for s,t in zip(sites,times):
                mask=((1<<(2*r["tau"][s]))-1)<<t
                if masks[-1]&mask:
                    break
                for h in reversed(range(1,M)):
                    masks[h]|=masks[h-1]&mask
                masks[0]|=mask
            else:
                sol=realize(r,assign,dict(zip(sites,times)))
                if sol is not None and sol["D"][-1]<=horizon:
                    key=tuple(masks)
                    cost=sol["D"][-1]-r["base"]
                    profiles[key]=min(profiles.get(key,cost),cost)
    return profiles

def oracle(I):
    if len(I["routes"])!=2:
        raise ValueError("Two-team exhaustive check")
    left,right=[local_profiles(r,I["Q"],I["H"],I["M"]) for r in I["routes"]]
    # Occupancy exceeds M iff left >= l and right >= M+1-l on a common tick.
    return min(a+b for ma,a in left.items() for mb,b in right.items()
               if not any(x&y for x,y in zip(ma,reversed(mb))))
