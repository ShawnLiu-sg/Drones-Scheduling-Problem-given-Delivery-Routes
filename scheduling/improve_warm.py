"""Fast feasible batch and dispatch-order search for a fixed FSTSP fleet."""
from bisect import insort
from collections import Counter
from dataclasses import dataclass
from math import inf
from pathlib import Path
from random import Random
import sys
from time import perf_counter

from scheduling.study import model, reference


def _close(dist):
    for k, row in enumerate(dist):
        for other in dist:
            lag = other[k]
            if lag == -inf:
                continue
            for j, tail in enumerate(row):
                if other[j] < lag+tail:
                    other[j] = lag+tail


@dataclass
class Pattern:
    assign: tuple
    sites: list
    early: list
    lags: list
    tail: list
    duration: list
    delay: int
    work: int


class Route:
    def __init__(self, rt, capacity):
        self.rt, self.capacity = rt, capacity
        self.allowed = [model.eligible(rt,p) for p in rt["parcels"]]
        n, m = len(rt["service"]), len(rt["sorties"])
        self.end = 2*n-1
        self.dist = [[0 if i==j else -inf for j in range(2*n+2*m)]
                     for i in range(2*n+2*m)]
        def arc(a,b,lag):
            self.dist[a][b] = max(self.dist[a][b],lag)
        for i in range(1,n):
            travel = rt["travel"][i-1]
            arc(n+i-1,i,travel)
            arc(i,n+i-1,-travel)
            arc(i,n+i,rt["service"][i])
        for j, f in enumerate(rt["sorties"]):
            launch, recovery = 2*n+j,2*n+m+j
            for evt,i in ((launch,f["l"]),(recovery,f["b"])):
                arc(i,evt,0)
                arc(evt,n+i,0)
            arc(launch,recovery,f["F"])
            arc(recovery,launch,-f["E"])
            if j:
                arc(recovery-1,launch,0)
        _close(self.dist)
        self.base = model.timetable([rt],[[]])[0]["A"]
        self.n = n
        self.cache = {}

    def pattern(self, assignment):
        if assignment in self.cache:
            return self.cache[assignment]
        rt, parcels = self.rt,self.rt["parcels"]
        sites=sorted(set(assignment))
        duration=[2*rt["tau"][i] for i in sites]
        early, ends, lags = [], [], []
        for i in sites:
            items = [p for p,loc in zip(parcels,assignment) if loc==i]
            if sum(p["w"] for p in items)>self.capacity:
                return None
            tau=rt["tau"][i]
            early.append(max(0,self.base[i]-tau,max(p["r"] for p in items)))
            service=rt["service"][i] if any(p["kind"]=="truck" and p["dest"]==i for p in items) else 0
            effects=[(self.n+i,tau+service)]
            effects.extend((2*self.n+p["dest"],tau) for p in items if p["kind"]=="sortie")
            ends.append(max(lag+self.dist[e][self.end] for e,lag in effects)-rt["base"])
            lags.append([max(lag+self.dist[e][j] for e,lag in effects)-rt["tau"][j] for j in sites])
        for j,row in enumerate(lags):
            row[j]=max(0,row[j])
        _close(lags)
        if any(row[j]>0 for j,row in enumerate(lags)):
            return None
        early = [max(t+row[j] for t,row in zip(early,lags)) for j in range(len(sites))]
        tail = [max(lag+end for lag,end in zip(row,ends)) for row in lags]
        delay=max([0]+[t+d for t,d in zip(early,tail)])
        pattern=Pattern(assignment,sites,early,lags,tail,duration,delay,sum(duration))
        self.cache[assignment]=pattern
        return pattern

    def neighbor(self, pattern, rng):
        if not pattern.assign:
            return pattern
        pos=rng.randrange(len(pattern.assign))
        assignment=list(pattern.assign)
        choices=[i for i in self.allowed[pos] if i != assignment[pos]]
        occupied=[i for i in choices if i in pattern.sites]
        if occupied and rng.random()<.66:
            choices=occupied
        if not choices:
            return pattern
        if rng.random()<.20:
            site=min(rng.choices(choices,k=3),key=lambda i:self.rt["tau"][i])
        else:
            site=rng.choice(choices)
        if rng.random()<.15:
            other=[j for j,place in enumerate(assignment) if place==site
                   and assignment[pos] in self.allowed[j]]
            if other:
                assignment[rng.choice(other)]=assignment[pos]
        assignment[pos]=site
        return self.pattern(tuple(assignment))


def _schedule(patterns, fleet, horizon, order, priorities=None):
    """Place each job on the first available machine, propagating route lags."""
    free=[0]*fleet
    pos=[0]*len(patterns)
    ready=[p.early.copy() for p in patterns]
    ends=[p.delay for p in patterns]
    starts=[[] for p in patterns]
    if order is None:
        order=[]
        count=sum(len(p.sites) for p in patterns)
    else:
        count=len(order)
    for step in range(count):
        avail=free[0]
        if priorities is None:
            k=order[step]
        else:
            w_slack,w_work,bias=priorities
            def rank(k):
                p=patterns[k]
                j=pos[k]
                t=max(avail,ready[k][j])
                slack=ends[k]-t-p.tail[j]
                return t+w_slack*slack+w_work*p.duration[j]+bias[k]
            k=min((k for k,p in enumerate(patterns) if pos[k]<len(p.sites)),key=rank)
            order.append(k)
        j=pos[k]
        pattern=patterns[k]
        t=max(avail,ready[k][j])
        finish=t+pattern.duration[j]
        if finish>horizon:
            return None
        for i,lag in enumerate(pattern.lags[j]):
            needed=t+lag
            if needed>ready[k][i]:
                if i<j and needed>starts[k][i]:
                    return None
                ready[k][i]=needed
        ends[k]=max(ends[k],t+pattern.tail[j])
        starts[k].append(t)
        pos[k]+=1
        free.pop(0)
        insort(free,finish)
    return sum(ends),order,starts


def improve(instance,time_limit=2.0,seed=13):
    tick=perf_counter()
    deadline=tick+time_limit
    rng=Random(seed)
    routes=instance["routes"]
    fleet=instance["M"]
    warm=instance.get("warm")
    best=warm
    objective=reference.validate(instance,warm) if warm is not None else None
    examined=0
    def result():
        return dict(objective=objective,routes=best,stats=dict(
            moves=examined,seconds=perf_counter()-tick))
    if not warm or objective==0 or time_limit<=0 or fleet<=0:
        return result()
    choices=[]
    pools=[]
    congestion=len(routes)/fleet
    for k,rt in enumerate(routes):
        if perf_counter()>=deadline:
            return result()
        explorer=Route(rt,instance["Q"])
        start=explorer.pattern(tuple(warm[k]["assign"]))
        choices.append(explorer)
        patterns=[start]
        for trial in range(100):
            if perf_counter()>=tick+time_limit*.32:
                break
            scale=rng.choice((.08,.30,.8,1.8,3.0))*congestion
            parent=min(rng.choices(patterns,k=min(12,len(patterns))),
                       key=lambda p:p.delay+scale*p.work)
            child=explorer.neighbor(parent,rng)
            if child is not None and child is not parent:
                patterns.append(child)
        pools.append(patterns)

    best_state=None
    current=None
    def consider(candidate, patterns):
        nonlocal objective,best,best_state
        if candidate is None:
            return
        value,order,st=candidate
        if value>=objective:
            return
        fixed=[dict(zip(p.sites,s)) for p,s in zip(patterns,st)]
        witness=model.timetable(routes,[p.assign for p in patterns],fixed=fixed)
        if witness is None or any(w["D"][-1]>instance["H"] for w in witness):
            return
        actual=reference.validate(instance,witness)
        if actual!=value:
            raise ArithmeticError(("Projected route delays differ from replay",value,actual))
        objective,best=actual,witness
        best_state=(patterns,order,value)
    initial=[p[0] for p in pools]
    baseline_order=[k for _,k,_ in sorted((t,k,int(i)) for k,w in enumerate(warm) for i,t in w["starts"].items())]
    init=_schedule(initial,fleet,instance["H"],baseline_order)
    if init:
        current=(initial,baseline_order,init[0])
        consider(init,initial)
    while perf_counter()<deadline-.003:
        examined+=1
        if current is None or examined<=16 or examined%120==0:
            scale=rng.choice((.08,.3,.6,1.1,2.0))*congestion
            pats=[min(rng.choices(pool,k=min(25,len(pool))),
                      key=lambda p:p.delay+scale*p.work) for pool in pools]
            bias=[rng.uniform(-3,3) for _ in routes]
            cand=_schedule(pats,fleet,instance["H"],None,
                           (rng.uniform(0,.6),rng.uniform(-.2,.7),bias))
        else:
            if examined%50==0 and best_state:
                current=best_state
            pats,order,v=current
            order=order.copy()
            if rng.random()<.42:
                k=rng.randrange(len(routes))
                old=pats[k]
                new=choices[k].neighbor(old,rng) if rng.random()<.72 else rng.choice(pools[k])
                if new is None or new is old:
                    continue
                pats=pats.copy()
                pats[k]=new
                delta=len(new.sites)-len(old.sites)
                for _ in range(-delta):
                    order.pop(rng.choice([i for i,team in enumerate(order) if team==k]))
                for _ in range(delta):
                    order.insert(rng.randrange(len(order)+1),k)
            elif len(order)>1:
                i=rng.randrange(len(order))
                if rng.random()<.80:
                    j=max(0,min(len(order)-1,i+rng.randint(-6,6)))
                else:
                    j=rng.randrange(len(order))
                order.insert(j,order.pop(i))
            cand=_schedule(pats,fleet,instance["H"],order)
        if cand is None:
            continue
        consider(cand,pats)
        if current is None or cand[0]<=current[2] or (
                cand[0]<=current[2]+rng.choice((1,2,3,5)) and rng.random()<.035):
            current=(pats,cand[1],cand[0])
    return result()
