"""Optional exact column restrictions and globally valid site-pair cuts for D."""
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from time import perf_counter

import gurobipy as gp


@dataclass
class PairCut:
    atoms: dict

    def hits(self, team, starts):
        return team in self.atoms and not self.atoms[team].isdisjoint(starts)


class Strengthening:
    def __init__(self, local=True, pairs=False, cut_limit=24):
        self.local = local
        self.pairs = pairs
        self.cut_limit = cut_limit
        self.cuts = []
        self.mask_keys = set()
        self.sites = {}
        self.clocks = {}
        self.stats = dict(local_rows=0, pair_cuts=0, separation_passes=0,
                          separation_seconds=0., prepare_seconds=0., max_violation=0.)

    @property
    def separates(self):
        return bool(self.clocks) and len(self.cuts) < self.cut_limit

    def prepare(self, instance, pricers, slots):
        tick = perf_counter()
        fleet = instance["M"]
        backgrounds = []
        for k, p in enumerate(pricers):
            times = defaultdict(list)
            for i,t in p.y:
                times[i].append(t)
            sites = sorted(times)
            tau = p.rt["tau"]
            if self.local:
                for i,j in zip(sites,sites[1:]):
                    if sum(p.rt["travel"][i:j]) < abs(tau[i]-tau[j]):
                        raise ValueError("Route does not satisfy ground/flight time dominance")
                use = {i:p.m.getVarByName(f"used[{i}]") for i in sites}
                start = {i:p.m.getVarByName(f"release[{i}]") for i in sites}
                # With an unknown subset of used sites, impose the rank lag
                # whenever all sites in a subsequence of length M+1 are used.
                if len(p.allowed) > fleet:
                    for group in combinations(sites,fleet+1):
                        i,j = group[0],group[-1]
                        relax = max(times[i])+2*tau[i]
                        inactive = gp.quicksum(1-use[s] for s in group)
                        p.m.addConstr(start[j] >= start[i]+2*tau[i]-relax*inactive)
                        self.stats["local_rows"] += 1
            if self.pairs:
                self.sites[k] = {i:tuple((i,t) for t in ts) for i,ts in times.items()}
                # A parcel with only one eligible site forces its use. For all
                # its feasible starts the intersection below is occupied.
                required = {choices[0] for choices in p.allowed if len(choices)==1}
                background = [0]*slots
                for i in required.intersection(sites):
                    for t in range(max(times[i]),min(times[i])+2*tau[i]):
                        background[t] += 1
                backgrounds.append(background)
        if self.pairs:
            totals = list(map(sum,zip(*backgrounds)))
            for a,b in combinations(range(len(pricers)),2):
                clock = [0]
                for t,total in enumerate(totals):
                    residual = fleet - (total-backgrounds[a][t]-backgrounds[b][t])
                    clock.append(clock[-1]+int(residual<=1))
                if clock[-1]:
                    self.clocks[a,b] = clock
        self.durations = [p.rt["tau"] for p in pricers]
        self.stats["prepare_seconds"] = perf_counter()-tick

    def separate(self, weights, deadline=None):
        tick = perf_counter()
        self.stats["separation_passes"] += 1
        mass = defaultdict(lambda:defaultdict(float))
        for team,starts,x in weights:
            for i,t in starts:
                mass[team,i][t] += x
        candidates = []
        for (a,b),clock in self.clocks.items():
            for i in self.sites[a]:
                da = 2*self.durations[a][i]
                left = [(clock[t],clock[t+da],w) for t,w in mass[a,i].items()
                        if clock[t]<clock[t+da]]
                if not left:
                    continue
                for j in self.sites[b]:
                    db = 2*self.durations[b][j]
                    right = [(clock[t],clock[t+db],w) for t,w in mass[b,j].items()
                             if clock[t]<clock[t+db]]
                    if not right:
                        continue
                    best = None
                    for u in sorted({s for s,e,w in left}):
                        if deadline is not None and perf_counter()>=deadline:
                            self.stats["separation_seconds"] += perf_counter()-tick
                            return 0
                        for v in sorted({e for s,e,w in left}):
                            val = sum(w for s,e,w in left if s<=u and e>=v)
                            val += sum(w for s,e,w in right if s<v and e>u)
                            if val>1+1e-7 and (best is None or val>best[0]):
                                best = (val,u,v)
                    if best:
                        candidates.append((*best,a,i,b,j))
        new = 0
        for val,u,v,a,i,b,j in sorted(candidates,reverse=True):
            clk = self.clocks[a,b]
            da,db = 2*self.durations[a][i],2*self.durations[b][j]
            la = frozenset((i,t) for _,t in self.sites[a][i]
                           if clk[t]<clk[t+da] and clk[t]<=u and clk[t+da]>=v)
            lb = frozenset((j,t) for _,t in self.sites[b][j]
                           if clk[t]<clk[t+db] and clk[t]<v and clk[t+db]>u)
            sig = a,la,b,lb
            if sig not in self.mask_keys:
                self.mask_keys.add(sig)
                self.cuts.append(PairCut({a:la,b:lb}))
                new += 1
                self.stats["max_violation"] = max(self.stats["max_violation"],val-1)
            if new>=4 or len(self.cuts)>=self.cut_limit:
                break
        self.stats["pair_cuts"] = len(self.cuts)
        self.stats["separation_seconds"] += perf_counter()-tick
        return new

    def price_terms(self, prices):
        terms = [defaultdict(float) for _ in self.sites]
        for cut,price in zip(self.cuts,prices):
            if price:
                for k,atoms in cut.atoms.items():
                    for atom in atoms:
                        terms[k][atom] += price
        return terms
