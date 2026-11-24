"""Fixed-batch timing and one-flight support ablations. Run as a module from the project root."""
import argparse
from bisect import bisect_right
from collections import deque
import csv
from itertools import accumulate
import json
from pathlib import Path
from random import Random
import sys
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB
from .. import d_price
from ..model import realize



NEG = float("-inf")


class Network:
    """Integral difference constraints for one fixed assignment and its hard fleet cap."""
    def __init__(self, route, assignment, fleet, horizon, flight_horizon=None):
        self.names = []
        self.ids = {}
        self.arcs = []

        def arc(u, v, lag):
            for name in (u, v):
                if name not in self.ids:
                    self.ids[name] = len(self.names)
                    self.names.append(name)
            self.arcs.append((self.ids[u], self.ids[v], lag))

        def eq(u, v, lag):
            arc(u, v, lag)
            arc(v, u, -lag)

        n = len(route["service"])
        eq("zero", ("A", 0), 0)
        eq("zero", ("D", 0), 0)
        for i in range(1, n):
            eq(("D", i-1), ("A", i), route["travel"][i-1])
            arc(("A", i), ("D", i), route["service"][i])
        for j, leg in enumerate(route["sorties"]):
            for kind, i in (("L", leg["l"]), ("U", leg["b"])):
                arc(("A", i), (kind, j), 0)
                arc((kind, j), ("D", i), 0)
            arc(("L", j), ("U", j), leg["F"])
            arc(("U", j), ("L", j), -leg["E"])
            if j:
                arc(("U", j-1), ("L", j), 0)
        self.stops = sorted(set(assignment))
        self.duration = [2*route["tau"][i] for i in self.stops]
        for i in self.stops:
            flight = route["tau"][i]
            if flight_horizon is not None:
                arc(("R", i), "zero", 2*flight - flight_horizon)
            arc(("A", i), ("R", i), -flight)
            arc(("R", i), ("D", i), flight)
        for job, i in zip(route["parcels"], assignment):
            arc("zero", ("R", i), job["r"])
            if job["kind"] == "sortie":
                arc(("R", i), ("L", job["dest"]), route["tau"][i])
            elif job["dest"] == i:
                arc(("R", i), ("D", i), route["tau"][i]+route["service"][i])
        for pos in range(len(self.stops)-fleet):
            arc(("R", self.stops[pos]), ("R", self.stops[pos+fleet]), self.duration[pos])
        arc(("D", n-1), "zero", -horizon)
        self.zero = self.ids["zero"]
        self.end = self.ids["D", n-1]
        self.starts = [self.ids["R", i] for i in self.stops]
        self.horizon = horizon

    def distances(self,allow_infeasible=False,deadline=None):
        n = len(self.names)
        d = [[NEG]*n for _ in range(n)]
        for i in range(n):
            d[i][i] = 0
        for i,j,w in self.arcs:
            d[i][j] = max(d[i][j],w)
        for k,row in enumerate(d):
            if deadline is not None and perf_counter()>=deadline:
                raise TimeoutError("difference closure")
            for r in d:
                w = r[k]
                if w != NEG:
                    for j,x in enumerate(row):
                        if w+x > r[j]:
                            r[j] = w+x
        if any(d[i][i] > 0 for i in range(n)):
            if allow_infeasible:
                return None
            raise ValueError("Infeasible base schedule in benchmark")
        return d


class Flow:
    def __init__(self,n):
        self.adj = [[] for _ in range(n)]

    def arc(self,u,v,cap):
        back,forward=len(self.adj[v]),len(self.adj[u])
        self.adj[u].append([v,back,cap])
        self.adj[v].append([u,forward,0])

    def cut(self,s,t,deadline=None):
        remaining_checks = 0
        def check():
            nonlocal remaining_checks
            if deadline is not None:
                if remaining_checks == 0:
                    if perf_counter() >= deadline:
                        raise TimeoutError("max-flow solve")
                    remaining_checks = 255
                else:
                    remaining_checks -= 1
        sys.setrecursionlimit(max(sys.getrecursionlimit(),min(20000,len(self.adj)+128)))
        while True:
            check()
            level = [-1]*len(self.adj)
            level[s] = 0
            q = deque([s])
            while q:
                check()
                u=q.popleft()
                for v,rev,cap in self.adj[u]:
                    if cap>1e-12 and level[v]<0:
                        level[v]=level[u]+1
                        q.append(v)
            if level[t]<0:
                return {v for v,lev in enumerate(level) if lev>=0}
            idx=[0]*len(self.adj)
            def push(u,limit):
                if u==t:
                    return limit
                while idx[u]<len(self.adj[u]):
                    check()
                    edge=self.adj[u][idx[u]]
                    v,rev,cap=edge
                    if cap>1e-12 and level[v]==level[u]+1:
                        sent=push(v,min(limit,cap))
                        if sent>1e-12:
                            edge[2]-=sent
                            self.adj[v][rev][2]+=sent
                            return sent
                    idx[u]+=1
                return 0
            while push(s,float("inf"))>1e-12:
                pass


def sparse_timing(net, dist, prices, terminal_scale=1, deadline=None):
    thresholds=[]
    allnodes=[]
    weights=[]
    edges=[]
    terminal={}
    c0=dist[net.zero][net.end]
    for i,ell in zip(net.starts,net.duration):
        if deadline is not None and perf_counter()>=deadline:
            raise TimeoutError("sparse threshold construction")
        lo,hi=dist[net.zero][i],-dist[i][net.zero]
        jumps={}
        for slot,p in prices.items():
            for a,delta in ((slot-ell+1,p),(slot+1,-p)):
                if lo<a<=hi:
                    jumps[a]=jumps.get(a,0)+delta
        times,ids=[],[]
        for a,dq in sorted(jumps.items()):
            if dq==0:
                continue
            v=len(weights)
            if ids:
                edges.append((v,ids[-1]))
            ids.append(v)
            times.append(a)
            weights.append(dq)
            allnodes.append((i,a))
            earliest_end=a+dist[i][net.end]
            if earliest_end>c0:
                terminal[earliest_end]=None
        thresholds.append((i,times,ids))
    prev=c0
    for c in sorted(terminal):
        v=len(weights)
        terminal[c]=v
        if prev>c0:
            edges.append((v,terminal[prev]))
        weights.append((c-prev)*terminal_scale)
        prev=c
    for u,(i,a) in enumerate(allnodes):
        if u%128==0 and deadline is not None and perf_counter()>=deadline:
            raise TimeoutError("sparse graph construction")
        for j,times,ids in thresholds:
            if i!=j:
                idx=bisect_right(times,a+dist[i][j])-1
                if idx>=0:
                    edges.append((u,ids[idx]))
        c=a+dist[i][net.end]
        if c>c0:
            edges.append((u,terminal[c]))
    src=len(weights)
    f=Flow(src+2)
    infinite=sum(abs(w) for w in weights)+1
    for u,v in edges:
        f.arc(u,v,infinite)
    for u,w in enumerate(weights):
        if w<0:
            f.arc(src,u,-w)
        if w>0:
            f.arc(u,src+1,w)
    selected=f.cut(src,src+1,deadline)
    times=dist[net.zero].copy()
    for v,(i,a) in enumerate(allnodes):
        if v in selected:
            for j,x in enumerate(dist[i]):
                times[j]=max(times[j],a+x)
    starts=[times[i] for i in net.starts]
    value=times[net.end]*terminal_scale+sum(p for s,e in zip(starts,net.duration)
                             for t,p in prices.items() if s<=t<s+e)
    if deadline is not None and perf_counter()>deadline:
        raise TimeoutError("sparse solution recovery")
    return value,starts,len(weights),len(edges)


class TimeMIP:
    """Time-indexed reference. Occupancy rows use original slot counts for support checks."""
    def __init__(self,net,dist,env,profile=False,fleet=None):
        self.model=m=gp.Model(env=env)
        m.Params.Threads=1
        m.Params.MIPGap=0
        m.Params.MIPGapAbs=0
        self.net=net
        clock=m.addVars(len(net.names),lb=[dist[net.zero][i] for i in range(len(net.names))],
                        ub=[-row[net.zero] for row in dist])
        for u,v,d in net.arcs:
            m.addConstr(clock[v]>=clock[u]+d)
        self.end=clock[net.end]
        self.starts=[]
        occup=[]
        if profile:
            occup=[[] for _ in range(net.horizon+max(net.duration))]
        for i,ell in zip(net.starts,net.duration):
            ts=range(dist[net.zero][i],1-dist[i][net.zero])
            y=m.addVars(ts,vtype=GRB.BINARY)
            m.addConstr(y.sum()==1)
            m.addConstr(clock[i]==gp.quicksum(t*y[t] for t in ts))
            self.starts.append(y)
            if profile:
                for t in ts:
                    for j in range(t,t+ell):
                        occup[j].append(y[t])
        self.rows=[m.addConstr(gp.quicksum(vs)<=fleet) for vs in occup]

    def timing(self,prices):
        prefix=list(accumulate((prices.get(t,0) for t in range(self.net.horizon+max(self.net.duration))),initial=0))
        self.model.setObjective(self.end+gp.quicksum((prefix[t+ell]-prefix[t])*v for y,ell in zip(self.starts,self.net.duration) for t,v in y.items()))
        self.model.optimize()
        if self.model.Status!=GRB.OPTIMAL:
            raise RuntimeError(self.model.Status)
        return self.model.ObjVal,[next(t for t,v in y.items() if v.X>.5) for y in self.starts]

    def support(self,times,ell,fleet):
        found=set()
        for s in times:
            if s==0:
                for j in range(ell):
                    self.rows[j].RHS=fleet-1
            else:
                self.rows[s-1].RHS=fleet
                self.rows[s+ell-1].RHS=fleet-1
            self.model.optimize()
            if self.model.Status==GRB.OPTIMAL:
                found.add(s)
            elif self.model.Status!=GRB.INFEASIBLE:
                raise RuntimeError(self.model.Status)
        return found


def append_arc(d,u,v,w):
    if w<=d[u][v]:
        return True
    if d[v][u]+w>0:
        return False
    row=d[v]
    for r in d:
        shift=r[u]+w
        if shift!=NEG:
            for j,value in enumerate(row):
                if shift+value>r[j]:
                    r[j]=shift+value
    return True


def project(net,dist,ell,fleet,max_start):
    # Project out all events except the starts and the origin before inserting the probe.
    ordered=[net.zero]+net.starts
    base=[[dist[i][j] for j in ordered]+[NEG] for i in ordered]
    base.append([NEG]*len(ordered)+[0])
    t=len(ordered)
    append_arc(base,0,t,0)
    append_arc(base,t,0,-max_start)
    b=len(net.starts)
    ranges=[]
    for u in range(b+1):
        for v in (range(u,u+1) if fleet==1 else range(u,b+1)):
            d=[row.copy() for row in base]
            if u and not append_arc(d,u,t,net.duration[u-1]):
                continue
            if v<b and not append_arc(d,t,v+1,ell):
                continue
            if all(append_arc(d,j+1,j+fleet,net.duration[j]) for j in range(u,v-fleet+1)):
                ranges.append((d[0][t],-d[t][0]))
    merged=[]
    for lo,hi in sorted(ranges):
        if merged and lo<=merged[-1][1]+1:
            merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
        else:
            merged.append((lo,hi))
    return merged
