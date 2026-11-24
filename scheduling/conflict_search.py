"""Exact conflict-based search (CBS) with fixed-route time-network projection."""
from bisect import bisect_left
from collections import Counter, defaultdict
from heapq import heappop, heappush
from itertools import count
from math import inf
from pathlib import Path
import sys
from time import perf_counter

from scheduling import liu_flow as liu
from scheduling.improve_warm import Route


class LocalOracle:
    """Integral parcel and departure search on one STN with site/time restrictions."""
    def __init__(self,route,capacity,horizon,deadline):
        self.rt = route
        self.capacity = capacity
        self.horizon = horizon
        self.deadline = deadline
        projection = Route(route,capacity)
        self.dist = projection.dist
        self.n = len(route["tau"])
        self.end = projection.end
        base = liu.schedules.timetable([route],[[]])[0]
        self.initial = base["A"]+base["D"]+base["L"]+base["U"]
        self.eligible = [[i for i in sites if 0<i<self.n-1] for sites in projection.allowed]
        self.calls=0
        self.states=0
        self.cache={}

    def _effects(self,p,i):
        parcel=self.rt["parcels"][p]
        tau=self.rt["tau"][i]
        effects=[(self.n+i,tau)]
        if parcel["kind"]=="truck" and parcel["dest"]==i:
            effects[0]=(self.n+i,tau+self.rt["service"][i])
        elif parcel["kind"]=="sortie":
            effects.append((2*self.n+parcel["dest"],tau))
        return effects

    def feasible(self,delay,forbidden,prefix):
        self.calls+=1
        rt=self.rt
        cap=rt["base"]+delay
        # Each fixed handoff and assignment adds lower bounds to site departure
        # and (if required) launch; receipt also upper-bounds truck arrival.
        floor=self.initial
        options=[]
        for p,parcel in enumerate(rt["parcels"]):
            if parcel["w"]>self.capacity:
                return None
            candidate=[]
            for i in self.eligible[p]:
                tau=rt["tau"][i]
                effects=self._effects(p,i)
                first=max(0,parcel["r"],floor[i]-tau)
                last=min(self.horizon-2*tau,
                         *(cap-lag-self.dist[node][self.end] for node,lag in effects))
                for t in range(first,int(last)+1):
                    if t % 256==0 and perf_counter()>=self.deadline:
                        raise TimeoutError
                    ticks=forbidden.get(i,())
                    idx=bisect_left(ticks,t)
                    if idx<len(ticks) and ticks[idx]<t+2*tau:
                        continue
                    conflict=prefix[t+2*tau]-prefix[t] if prefix is not None else 0
                    candidate.append((conflict,t+tau,i,t,effects))
            if not candidate:
                return None
            candidate.sort()
            options.append(candidate)
        ordering=sorted(range(len(options)), key=lambda p:(len(options[p]),-rt["parcels"][p]["w"]))
        assign=[None]*len(options)
        starts={}
        loads=Counter()

        def search(depth,earliest):
            self.states+=1
            if self.states%64==0 and perf_counter()>=self.deadline:
                raise TimeoutError
            if depth==len(ordering):
                return liu.schedules.realize(rt,assign,starts)
            p=ordering[depth]
            weight=rt["parcels"][p]["w"]
            # Reusing an existing dispatch consumes no additional fleet time.
            for existing in (True,False):
                for _,receipt,i,t,effects in options[p]:
                    if (i in starts)!=existing:
                        continue
                    if existing and starts[i]!=t or loads[i]+weight>self.capacity:
                        continue
                    if earliest[i]>receipt:
                        continue
                    after=earliest
                    for node,lag in effects:
                        value=t+lag
                        if value>after[node]:
                            after=[max(v,value+w) for v,w in zip(after,self.dist[node])]
                    if after[self.end]>cap or after[0]>0 or after[self.n]>0 or after[i]>receipt:
                        continue
                    if any(after[j]>s+rt["tau"][j] for j,s in starts.items()):
                        continue
                    assign[p]=i
                    starts[i]=t
                    loads[i]+=weight
                    solution=search(depth+1,after)
                    if solution is not None:
                        return solution
                    loads[i]-=weight
                    if not existing:
                        del starts[i]
            return None
        return search(0,floor)

    def minimum(self,key,lower,upper,prefix):
        known,plan=self.cache.get(key,(lower,None))
        if plan is not None:
            return "optimal",known,plan
        lower=max(lower,known)
        forbidden=defaultdict(list)
        for i,t in sorted(key):
            forbidden[i].append(t)
        while lower<=upper:
            if perf_counter()>=self.deadline:
                return "limit",lower,None
            try:
                plan=self.feasible(lower,forbidden,prefix)
            except TimeoutError:
                return "limit",lower,None
            if plan is not None:
                value=plan["D"][-1]-self.rt["base"]
                if value!=lower:
                    raise ArithmeticError(("Local search skipped a feasible delay",lower,value))
                self.cache[key]=(value,plan)
                return "optimal",value,plan
            lower+=1
            self.cache[key]=(lower,None)
        return "excluded",lower,None


def _conflict(routes,plans,fleet):
    sweep=defaultdict(lambda:[[],[]])
    for k,(rt,plan) in enumerate(zip(routes,plans)):
        for site,start in plan["starts"].items():
            site=int(site)
            sweep[start][1].append((k,site))
            sweep[start+2*rt["tau"][site]][0].append((k,site))
    active=set()
    for t,(end,start) in sorted(sweep.items()):
        active.difference_update(end)
        active.update(start)
        if len(active)>fleet:
            return t,sorted(active)[:fleet+1]
    return None


def _congestion(instance,plans,excluded):
    if not any(p is not None for k,p in enumerate(plans) if k!=excluded):
        return None
    horizon=instance["H"]
    events=Counter()
    for k,(rt,p) in enumerate(zip(instance["routes"],plans)):
        if p is None or k==excluded:
            continue
        for site,t in p["starts"].items():
            events[t]+=1
            events[t+2*rt["tau"][int(site)]]-=1
    total=0
    prefix=[0]
    for t in range(horizon):
        total+=events[t]
        prefix.append(prefix[-1]+int(total>=instance["M"]))
    return prefix


def solve(instance,time_limit=120,threads=1,*,progress=None,oracle_class=LocalOracle):
    """Solve exactly on the integer time grid; return frontier LB at timeout."""
    tick=perf_counter()
    deadline=tick+time_limit
    if "H" not in instance:
        return liu.solve(instance,time_limit,threads,progress=progress)
    routes=instance["routes"]
    best=instance.get("warm")
    upper=liu.validate(instance,best) if best is not None else inf
    lower=0
    oracles=[]
    nodes=0
    branches=0
    def publish():
        if progress:
            progress(upper if upper<inf else None,lower)
    def result(status):
        publish()
        return dict(status=status,objective=upper if upper<inf else None,
                    bound=None if status=="infeasible" else lower,routes=best,origin="algorithm",
                    stats=dict(nodes=nodes,branches=branches,local_decisions=sum(x.calls for x in oracles),
                               local_states=sum(x.states for x in oracles),seconds=perf_counter()-tick))
    baselines=liu.schedules.timetable(routes,[[] for _ in routes])
    if baselines is None:
        return result("infeasible")
    lbs=[]
    for rt,base in zip(routes,baselines):
        latest=[min(liu._latest(rt,p),len(rt["tau"])-2) for p in rt["parcels"]]
        v=liu._team_lower_bound(rt,base,liu._road_tail(rt),latest,instance["Q"])
        if v is None or rt["base"]+v>instance["H"]:
            return result("infeasible")
        lbs.append(v)
    lower=sum(lbs)
    if lower==upper:
        return result("optimal")
    publish()
    for rt in routes:
        if perf_counter()>=deadline:
            return result("limit")
        oracles.append(oracle_class(rt,instance["Q"],instance["H"],deadline))
    seq=count()
    empty=tuple(frozenset() for _ in routes)
    frontier=[(lower,0,next(seq),empty,lbs,[None]*len(routes),list(range(len(routes))))]
    visited={empty}

    while frontier:
        bound,negdepth,_,prohibitions,inherited,proposals,dirty=heappop(frontier)
        lower=min(upper,bound)
        if lower==upper:
            return result("optimal")
        if perf_counter()>=deadline:
            return result("limit")
        nodes+=1
        local=inherited.copy()
        plans=proposals.copy()
        excluded=False
        for k in dirty:
            cost_cap=instance["H"]-routes[k]["base"]
            if upper<inf:
                cost_cap=min(cost_cap,int(upper)-1-sum(local)+local[k])
            status,v,plan=oracles[k].minimum(
                prohibitions[k],local[k],cost_cap,_congestion(instance,plans,k))
            local[k]=v
            if status=="limit":
                lower=min(upper,sum(local),frontier[0][0] if frontier else inf)
                return result("optimal" if lower==upper else "limit")
            if status=="excluded":
                excluded=True
                break
            plans[k]=plan
        if excluded or sum(local)>=upper:
            continue
        value=sum(local)
        if frontier and value>frontier[0][0]:
            heappush(frontier,(value,negdepth,next(seq),prohibitions,local,plans,[]))
            continue
        lower=value
        publish()
        collision=_conflict(routes,plans,instance["M"])
        if collision is None:
            actual=liu.validate(instance,plans)
            if actual!=value:
                raise ArithmeticError(("Local and joint cost differ",value,actual))
            best,upper=plans,actual
            return result("optimal")
        t,flights=collision
        for k,i in flights:
            forbidden=list(prohibitions)
            forbidden[k]=forbidden[k] | {(i,t)}
            key=tuple(forbidden)
            if key in visited:
                continue
            visited.add(key)
            branches+=1
            heappush(frontier,(value,negdepth-1,next(seq),key,local,plans,[k]))
    if best is not None:
        lower=upper
        return result("optimal")
    return result("infeasible")
