"""End-to-end pricing ablation with all feasible batch assignments and exact integer cuts."""
import argparse
import csv
from itertools import product
import json
import math
from pathlib import Path
from random import Random
from time import perf_counter

import gurobipy as gp
from .. import d_price,reference

from ..model import eligible,realize
from .temporal import Network,sparse_timing


def integer_prices(prefix):
    # Work from the exact binary rational representation of the solver's prefix array.
    fractions=[v.as_integer_ratio() for v in prefix]
    scale=max(b for a,b in fractions)
    sums=[a*(scale//b) for a,b in fractions]
    return scale,{t:b-a for t,(a,b) in enumerate(zip(sums,sums[1:])) if b>a}


class EnumerationPricing:
    def __init__(self,capacity,fleet=None):
        self.fleet=fleet
        self.capacity=capacity
        self.cache={}
        self.enumerated=0
        self.feasible=0
        self.oracle_calls=0
        self.complete_calls=0
        self.fallback_calls=0
        self.fixed_solves=0
        self.oracle_seconds=0
        self.build_seconds=0

    def patterns(self,pricer,deadline):
        if pricer in self.cache:
            return self.cache[pricer]
        tick=perf_counter()
        rt=pricer.rt
        choices=[eligible(rt,p) for p in rt["parcels"]]
        complete=[]
        for assignment in product(*choices):
            if perf_counter()>=deadline:
                self.build_seconds+=perf_counter()-tick
                return None
            self.enumerated+=1
            load={}
            for stop,parcel in zip(assignment,rt["parcels"]):
                load[stop]=load.get(stop,0)+parcel["w"]
            if any(w>self.capacity for w in load.values()):
                continue
            fleet = self.fleet if self.fleet is not None else len(assignment)+1
            net=Network(rt,assignment,fleet,pricer.horizon,pricer.flight_horizon)
            try:
                dist=net.distances(allow_infeasible=True,deadline=deadline)
            except TimeoutError:
                self.build_seconds+=perf_counter()-tick
                return None
            if dist is not None:
                complete.append((assignment,net,dist))
                self.feasible+=1
        self.cache[pricer]=complete
        self.build_seconds+=perf_counter()-tick
        return complete

    def price(self,pricer,dual,prefix,deadline):
        t0=perf_counter()
        self.oracle_calls+=1
        patterns=self.patterns(pricer,deadline)
        if patterns is None:
            self.oracle_seconds+=perf_counter()-t0
            return None,-math.inf
        scale,prices=integer_prices(prefix)
        best=None
        witness=None
        completed=True
        for assignment,net,distance in patterns:
            if perf_counter()>=deadline:
                completed=False
                break
            try:
                value,times,_,_=sparse_timing(net,distance,prices,terminal_scale=scale,deadline=deadline)
            except TimeoutError:
                completed=False
                break
            self.fixed_solves+=1
            if best is None or value<best:
                best=value
                witness=(assignment,dict(zip(net.stops,times)))
        self.complete_calls+=int(completed)
        if witness is None:
            answer=(None,math.inf if completed else -math.inf)
        else:
            assignment,starts=witness
            schedule=realize(pricer.rt,assignment,starts)
            if schedule is None or schedule["D"][-1]>pricer.horizon:
                raise ArithmeticError("Mincut timing produced a nonrealizable column")
            nd,dd=dual.as_integer_ratio()
            rc=((best-scale*pricer.base)*dd-scale*nd)/(scale*dd)
            rc_wit=(schedule["D"][-1]-pricer.base-dual+
                    sum(prefix[t+2*pricer.rt["tau"][i]]-prefix[t] for i,t in starts.items()))
            if abs(rc_wit-rc)>1e-6:
                raise ArithmeticError(("Wrong reduced cost",rc_wit,rc))
            # Sum of floating addition errors in building the master's prefix array,
            # uniformly for every schedule and every possible occupied interval.
            margin=len(prefix)*len(pricer.rt["parcels"])*math.ulp(prefix[-1])
            lb=math.nextafter(rc-margin,-math.inf) if completed else -math.inf
            answer=(schedule,rc_wit),lb
        self.oracle_seconds+=perf_counter()-t0
        return answer


def run(I,enabled,limit):
    engine=EnumerationPricing(I["Q"])
    old=d_price._Price.price
    def routed(self,dual,prefix,phase,fixed,deadline,cut_costs=None):
        if enabled and phase==2 and not fixed and not cut_costs:
            return engine.price(self,dual,prefix,deadline)
        engine.fallback_calls+=1
        return old(self,dual,prefix,phase,fixed,deadline,cut_costs)
    d_price._Price.price=routed
    try:
        start=perf_counter()
        res=d_price.solve(I,limit,1)
        elapsed=perf_counter()-start
    finally:
        d_price._Price.price=old
    stats=res["stats"]
    fields="enumerated feasible oracle_calls complete_calls fallback_calls fixed_solves oracle_seconds build_seconds".split()
    out={f:getattr(engine,f) for f in fields}
    out.update(method="enumerate_cut" if enabled else "original",status=res["status"],
               ub=res["objective"],lb=res["bound"],elapsed=elapsed,
               branches=stats["branches"],nodes=stats["nodes"],columns=stats["columns"],
               pricing_calls=stats["pricing_calls"],pricing_seconds=stats["pricing_seconds"])
    assert reference.validate(I,res["routes"])==res["objective"]
    return out,res["routes"]
