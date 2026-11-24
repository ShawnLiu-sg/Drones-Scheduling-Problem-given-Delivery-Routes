"""Capacity-conflict search with a fixed-assignment, earliest-dispatch oracle.

Every (site,tick) prohibition removes a closed interval of integer start times.
On a fixed batch assignment the remaining constraints are differences of times.
A monotone propagation jumps each dispatch to its next admissible start. It
therefore optimizes timing without enumerating combinations of dispatch times.
"""
from bisect import bisect_right
from collections import deque
from math import inf
from time import perf_counter

from scheduling import conflict_search as cbs


class Domain:
    def __init__(self, ticks, duration):
        self.left = []
        self.right = []
        for t in ticks:
            lo,hi = t-duration+1,t
            if self.right and lo <= self.right[-1]+1:
                self.right[-1] = hi
            else:
                self.left.append(lo)
                self.right.append(hi)

    def next(self, earliest):
        j = bisect_right(self.left,earliest)-1
        if j>=0 and earliest<=self.right[j]:
            return self.right[j]+1
        return earliest


class BatchOracle(cbs.LocalOracle):
    def __init__(self, route, capacity, horizon, deadline):
        super().__init__(route,capacity,horizon,deadline)
        # For each parcel/site choice, eliminate the static truck and sortie
        # clocks once. row[j] sends the source dispatch's lag to site j.
        self.options=[]
        self._previous=None
        self._domains=None
        for p,sites in enumerate(self.eligible):
            choices=[]
            for i in sites:
                tau=route["tau"][i]
                effects=self._effects(p,i)
                row=tuple(max(offset+self.dist[e][j] for e,offset in effects)-route["tau"][j]
                          for j in range(self.n))
                tail=max(offset+self.dist[e][self.end] for e,offset in effects)-route["base"]
                anchor=max(offset+max(self.dist[e][0],self.dist[e][self.n]) for e,offset in effects)
                first=max(0,route["parcels"][p]["r"],self.initial[i]-tau)
                last=min(horizon-2*tau,-anchor)
                choices.append((i,first,last,tail,row))
            self.options.append(choices)

    @staticmethod
    def _push(groups, clocks, index, minimum, domains):
        """Maintain the least dispatch vector after one new parcel decision."""
        times=clocks.copy()
        # The newly activated site also receives lags from every existing site.
        first=max([minimum, times[index]]
                  +[times[j]+group[3][index] for j,group in groups.items() if j!=index])
        times[index]=domains[index].next(first)
        if times[index]>groups[index][1]:
            return None
        queue=deque([index])
        while queue:
            j=queue.popleft()
            row=groups[j][3]
            for dest,group in groups.items():
                implied=times[j]+row[dest]
                if implied>times[dest]:
                    t=domains[dest].next(implied)
                    if t>group[1]:
                        return None
                    times[dest]=t
                    queue.append(dest)
        return times

    def feasible(self,delay,forbidden,prefix):
        self.calls+=1
        if forbidden is not self._previous:
            self._domains=[Domain(forbidden.get(i,()),2*tau) for i,tau in enumerate(self.rt["tau"])]
            self._previous=forbidden
        domains=self._domains
        rt=self.rt
        candidates=[]
        for p,opts in enumerate(self.options):
            if rt["parcels"][p]["w"]>self.capacity:
                return None
            possible=[]
            for i,first,last,tail,row in opts:
                first=domains[i].next(first)
                cap=min(last,delay-tail)
                if first<=cap:
                    occupancy=(prefix[first+2*rt["tau"][i]]-prefix[first]) if prefix is not None else 0
                    possible.append((occupancy,first+tail,rt["tau"][i],i,first,cap,tail,row))
            if not possible:
                return None
            candidates.append(sorted(possible))
        ordering=sorted(range(len(candidates)),key=lambda p:(len(candidates[p]),-rt["parcels"][p]["w"]))
        assign=[None]*len(candidates)
        groups={}
        def search(depth,clocks):
            self.states+=1
            if self.states%32==0 and perf_counter()>=self.deadline:
                raise TimeoutError
            if depth==len(ordering):
                starts={i:clocks[i] for i in groups}
                plan=cbs.liu.schedules.realize(rt,assign,starts)
                if plan is None:
                    raise ArithmeticError("Projected timing is inconsistent with the event graph")
                return plan
            p=ordering[depth]
            weight=rt["parcels"][p]["w"]
            for reused in (True,False):
                for _,_,_,i,first,cap,tail,row in candidates[p]:
                    old=groups.get(i)
                    if (old is not None)!=reused:
                        continue
                    if old is None:
                        merged=(weight,cap,tail,row)
                    else:
                        old_weight,old_cap,old_tail,old_row=old
                        if old_weight+weight>self.capacity:
                            continue
                        merged=(old_weight+weight,min(cap,old_cap),max(tail,old_tail),
                                tuple(map(max,row,old_row)))
                    groups[i]=merged
                    after=self._push(groups,clocks,i,first,domains)
                    if after is not None:
                        assign[p]=i
                        plan=search(depth+1,after)
                        if plan is not None:
                            return plan
                    if old is None:
                        del groups[i]
                    else:
                        groups[i]=old
            return None
        return search(0,[0]*self.n)


def solve(instance,time_limit=120,threads=1,*,progress=None):
    return cbs.solve(instance,time_limit,threads,progress=progress,oracle_class=BatchOracle)
