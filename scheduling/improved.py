from __future__ import annotations
from math import prod
from time import perf_counter
from collections import Counter

from scheduling.study.model import eligible
from scheduling.study.validation_d.temporal_global import EnumerationPricing
from scheduling.study.validation_d.temporal import Network
from itertools import product
from scheduling.domain_filter import filtered_choices
from scheduling.fast_network import fast_distances
from scheduling.projected_core import ProjectedRoute
from scheduling import engine


def capacity_assignments(route, choices, capacity, deadline):
    loads = [0]*len(route["tau"])
    assignment = [0]*len(choices)
    parcels = route["parcels"]
    def visit(p):
        if perf_counter() >= deadline:
            yield None
            return
        if p == len(choices):
            yield tuple(assignment)
            return
        w=parcels[p]["w"]
        for i in choices[p]:
            if loads[i]+w > capacity:
                continue
            loads[i]+=w
            assignment[p]=i
            yield from visit(p+1)
            loads[i]-=w
    yield from visit(0)

def timing_signature(route, assignment):
    # At any earlier station, receipt precedes departure and therefore every later launch.
    max_release={}
    same_truck=set()
    same_sortie=set()
    for p,i in zip(route["parcels"],assignment):
        max_release[i]=max(max_release.get(i,0),p["r"])
        if p["kind"]=="truck" and i==p["dest"]:
            same_truck.add(i)
        if p["kind"]=="sortie" and i==route["sorties"][p["dest"]]["l"]:
            same_sortie.add(p["dest"])
    return tuple(sorted(max_release.items())),tuple(sorted(same_truck)),tuple(sorted(same_sortie))

class ImprovedPricing(EnumerationPricing):
    def __init__(self,capacity,fleet,mode="all",threshold=2048):
        super().__init__(capacity,fleet)
        self.mode=mode
        self.threshold=threshold
        self.domains={}
        self.counts=Counter()
        self.cores={}
    def choices(self, pricer,deadline=None):
        if pricer not in self.domains:
            rt=pricer.rt
            if self.mode in ("filter","all","core","adaptive"):
                choices=filtered_choices(rt,pricer.horizon,pricer.flight_horizon,deadline)
            else:
                choices=[eligible(rt,p) for p in rt["parcels"]]
            self.domains[pricer]=choices
        return self.domains[pricer]
    def price(self,pricer,dual,prefix,deadline):
        try:
            domains=self.choices(pricer,deadline)
        except TimeoutError:
            return None,-float("inf")
        if self.mode=="adaptive" and prod(map(len,domains))>self.threshold:
            self.counts["switch_milp"]+=1
            return pricer.price(dual,prefix,2,{},deadline,None)
        return super().price(pricer,dual,prefix,deadline)
    def patterns(self,pricer,deadline):
        if pricer in self.cache: return self.cache[pricer]
        begin=perf_counter()
        rt=pricer.rt
        choices=self.choices(pricer)
        compressed=self.mode in ("signature","all","core","adaptive")
        projection=self.mode in ("closure","all","adaptive")
        self.counts["domain_product"]+=prod(map(len,choices))
        seen=set()
        complete=[]
        stream=capacity_assignments(rt,choices,self.capacity,deadline) if compressed else product(*choices)
        for row in stream:
            if perf_counter()>=deadline:
                self.build_seconds+=perf_counter()-begin
                return None
            self.enumerated+=1
            if compressed:
                key=timing_signature(rt,row)
                if key in seen:
                    self.counts["duplicates"]+=1
                    continue
                seen.add(key)
            else:
                load=Counter()
                for i,p in zip(row,rt["parcels"]):load[i]+=p["w"]
                if any(w>self.capacity for w in load.values()):continue
            self.counts["closures"]+=1
            fleet=self.fleet if self.fleet is not None else len(row)+1
            try:
                if self.mode in ("core","adaptive"):
                    if pricer not in self.cores:
                        self.cores[pricer]=ProjectedRoute(rt,pricer.horizon,pricer.flight_horizon,deadline)
                    net,dist=self.cores[pricer].pattern(row,fleet,deadline)
                else:
                    net=Network(rt,row,fleet,pricer.horizon,pricer.flight_horizon)
                    dist=(fast_distances(net,allow_infeasible=True,deadline=deadline) if projection else
                          net.distances(allow_infeasible=True,deadline=deadline))
            except TimeoutError:
                self.build_seconds+=perf_counter()-begin
                return None
            if dist is not None:
                complete.append((row,net,dist));self.feasible+=1
        self.build_seconds+=perf_counter()-begin
        self.cache[pricer]=complete
        return complete

_original=engine.Routing
def solve(I,mode,seconds,threshold=2048):
    if mode=="original":
        return engine.run(I,"full",seconds)
    class Router(_original):
        def __init__(self,I,hybrid,local):
            super().__init__(I,hybrid,local)
            self.enumerator=ImprovedPricing(I["Q"],I["M"] if local else None,mode,threshold)
        def stats(self):
            return dict(super().stats(),efficiency_counts=dict(self.enumerator.counts))
    engine.Routing=Router
    try:
        return engine.run(I,"full",seconds)
    finally:
        engine.Routing=_original
