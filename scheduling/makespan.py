"""Optimize the last return over the truck and shared replenishment fleets."""
from collections import defaultdict
from math import ceil,isfinite
from pathlib import Path
import re
import subprocess
from time import perf_counter
from scheduling.relaxation import _route_bound,warm_plan
from scheduling.study import model,reference
from scheduling import interval_exact as cp,certificate as cert,compact_flow as compact
from scheduling.endurance import read_cp_plan,digest
import gurobipy as gp
from gurobipy import GRB

def global_horizon(routes):
    tau = max(max(r["tau"]) for r in routes)
    count = 1 + sum(2*len(r["service"]) + 2*len(r["sorties"]) + len(r["parcels"]) for r in routes)
    lag = max([1, 2*tau] + [v for r in routes for v in
        [*r["travel"], max(r["service"])+max(r["tau"]),
         *(p["r"] for p in r["parcels"]), *(s["F"] for s in r["sorties"])]])
    return count*lag + 2*tau


def replay(inst, plan):
    truck = max(p["D"][-1] for p in plan)
    returns = [int(t)+2*r["tau"][int(i)] for r, p in zip(inst["routes"], plan, strict=True)
               for i,t in p["starts"].items()]
    last = max(returns, default=0)
    value = max(truck, last)
    reference.validate(dict(inst,H=value), plan)
    return int(value),int(truck),("both" if returns and truck==last else "truck" if truck>=last else "replenishment")


def energy_check(inst, H, raw):
    profiles = []
    for r, (options, _) in zip(inst["routes"], raw):
        windows = [{i:(f,l,min(s,H-2*r["tau"][i])) for i,(f,l,s) in ps.items()} for ps in options]
        profile = cert._classify(r, windows, H-r["base"], inst["Q"])
        if profile is None:
            return dict(kind="empty_parcel_or_forced_stop_domain")
        profiles.append(profile)
    points = {0,H}
    for forced, pieces in profiles:
        for l,e,p,*_ in [*forced,*pieces]:
            points.update(x for x in (l,e,l+p,e-p) if 0<=x<=H)
    ticks = sorted(points)
    if len(ticks)>25:
        ticks = [ticks[j*(len(ticks)-1)//24] for j in range(25)]
    windows = [(a,b) for j,a in enumerate(ticks) for b in ticks[j+1:]]
    demand = [0]*len(windows)
    for profile in profiles:
        for j, val in enumerate(cert._energy(profile, windows, inst["Q"])):
            demand[j] += val
    for used,(a,b) in zip(demand, windows):
        if used > inst["M"]*(b-a):
            return dict(kind="interval_energy",a=a,b=b,required=used,capacity=inst["M"]*(b-a))
    return None


def cp_trial(inst,H,lower,folder,mode,deadline,hint=None):
    begin = perf_counter()
    item = dict(method="interval_"+mode,threshold=H)
    local = dict(inst,H=H)
    program = cp._build(local,hint,None)
    if program is None:
        return dict(item,status="empty",reason="empty_interval_domain"),None,H+1
    source, stops = program
    lines = source.splitlines()
    assert len([x for x in lines if x.startswith("minimize(")])==1
    lines = [x for x in lines if not x.startswith("minimize(")]
    if mode=="optimize":
        ends = [f"D{k}_{len(r['tau'])-1}" for k,r in enumerate(inst["routes"])]
        ends += [f"endOf(j{k}_{i},0)" for k,row in enumerate(stops) for i in row]
        objective = [f"joint = intVar({lower}..{H});"]
        objective += [f"joint >= {x};" for x in ends]
        objective += ["minimize(joint);"]
        place = lines.index("startingPoint {") if "startingPoint {" in lines else len(lines)
        lines[place:place]=objective
    folder.mkdir(exist_ok=True,parents=True)
    path=folder/"model.cpo"
    path.write_text("\n".join(lines)+"\n",encoding="utf8")
    item.update(model=str(Path(folder.name)/path.name),model_sha256=digest(path))
    left = deadline-perf_counter()
    if left<.12:
        return dict(item,status="limit_build"),None,lower
    cmd=[str(cp.CP),"-c","read model.cpo","set Workers 1","set RandomSeed 13",
         f"set TimeLimit {max(.01,left-min(.7,.20*left)):.4f}",
         "set OptimalityTolerance 0","set RelativeOptimalityTolerance 0","optimize","display solution"]
    try:
        proc=subprocess.run(cmd,cwd=folder,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=left)
        log=proc.stdout.decode("utf8",errors="replace")
    except subprocess.TimeoutExpired as ex:
        log=(ex.stdout or b"").decode("utf8",errors="replace")
        (folder/"solve.log").write_text(log,encoding="utf8")
        return dict(item,status="limit_process"),None,lower
    plan=read_cp_plan(log,local,stops)
    if plan is not None:
        val,_,_=replay(local,plan)
        assert val<=H
    empty=bool(re.search(r"(?mi)^\s*(?:!\s*)?(?:Search completed, model (?:is infeasible|has no solution)|Model is infeasible)",log))
    bound=lower
    if empty:
        assert plan is None
        bound=H+1
    elif mode=="optimize":
        m=re.search(r"Best bound\s*:\s*(-?\d+(?:\.\d+)?)",log)
        if m:
            bound=max(lower,min(H+1,ceil(float(m.group(1))-1e-7)))
    if plan is not None:
        assert bound<=val
    (folder/"solve.log").write_text(log,encoding="utf8")
    if proc.returncode or re.search(r"(?mi)^.*(?:ILOG concert exception|Modeling error|Syntax error at)",log):
        raise RuntimeError(log[-2000:])
    return dict(item,status="empty" if empty else "witness" if plan else "limit",elapsed=perf_counter()-begin),plan,bound


def compact_trial(inst,H,base,folder,deadline):
    local=dict(inst,H=H)
    rows=[]
    estimate=0
    for r,b in zip(inst["routes"],base):
        lo,hi,_=compact._clock_graph(r,b,H-r["base"],H)
        if any(lo[x]>hi[x] for x in lo):
            return dict(method="compact",status="empty",threshold=H),None
        estimate+=sum(hi[x]-lo[x] for x in lo)
        if estimate>5000:
            return dict(method="compact",status="wide_use_interval_oracle",threshold=H),None
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag",0)
        env.start()
        with gp.Model("absolute_threshold",env=env) as m:
            m.Params.Threads=1
            m.Params.Seed=13
            m.Params.IntFeasTol=1e-8
            m.Params.FeasibilityTol=1e-8
            for r,b in zip(inst["routes"],base):
                latest=[min(compact.liu._latest(r,p),len(r["tau"])-2) for p in r["parcels"]]
                row=compact._team(m,local,r,b,latest,H-r["base"],0,deadline)
                if row is None:
                    return dict(method="compact",status="empty",threshold=H),None
                rows.append(row)
            ticks=defaultdict(gp.LinExpr)
            for stops,_,_,_ in rows:
                for dispatch in stops.values():
                    for t in range(dispatch.first,dispatch.last+2*dispatch.tau):
                        ticks[t]+=dispatch.at(t-2*dispatch.tau+1)-dispatch.at(t+1)
                if perf_counter()>=deadline:
                    raise TimeoutError
            for value in ticks.values():
                m.addLConstr(value<=local["M"])
            m.setObjective(0)
            remain=deadline-perf_counter()
            if remain<=.15:
                raise TimeoutError
            folder.mkdir(parents=True,exist_ok=True)
            m.Params.OutputFlag=1
            m.Params.LogToConsole=0
            m.Params.LogFile=str(folder/"gurobi.log")
            m.Params.TimeLimit=max(.01,remain-.15)
            m.optimize()
            info=dict(method="compact",threshold=H,status="empty" if m.Status==GRB.INFEASIBLE else "witness" if m.SolCount else "limit",
                      solver_status=m.Status,variables=m.NumVars,rows=m.NumConstrs,nodes=m.NodeCount,
                      log=str(Path(folder.name)/"gurobi.log"))
            if not m.SolCount:
                return info,None
            fixed,assign=[],[]
            for stops,loads,_,_ in rows:
                fixed.append({i:h.first+sum(v.X>.5 for v in h.tail)-1 for i,h in stops.items() if h.used.X>.5})
                assign.append([next(i for i,v in opts if v.X>.5) for opts in loads])
            plan=model.timetable(inst["routes"],assign,fixed=fixed)
            assert plan is not None and replay(local,plan)[0]<=H
            return info,plan


def result_status(lb,ub,global_empty=False):
    if global_empty:
        return "infeasible"
    if ub is None:
        return "limit_unknown"
    return "optimal" if lb==ub else "limit_feasible"


def solve(inst, cap, folder, initial=None):
    started=perf_counter()
    deadline=started+cap
    rt=inst["routes"]
    safe=min(global_horizon(rt), inst["H"])
    base=model.timetable(rt,[[] for _ in rt])
    local=[_route_bound(r) for r in rt]
    globally_empty=base is None or any(not isfinite(v) for v in local)
    lb=max((r["base"] for r in rt),default=0)
    ub=plan=val_time=truck=last=source=None
    history=[]
    if not globally_empty:
        lb=max(r["base"]+int(v) for r,v in zip(rt,local))
        for r in rt:
            for p in r["parcels"]:
                opts=model.eligible(r,p)
                if not opts or p["w"]>inst["Q"]:
                    globally_empty=True
                    break
                lb=max(lb,p["r"]+2*min(r["tau"][i] for i in opts))
    history.append(dict(method="individual_release_and_earliest_events",lower=lb,route_delay=local if not globally_empty else None))
    def save(candidate,why):
        nonlocal ub,plan,val_time,truck,last,source
        if candidate is not None:
            value,t,k=replay(inst,candidate)
            if value>safe:
                return
            assert value>=lb
            if perf_counter()<=deadline and (ub is None or value<ub):
                ub,truck,last,plan,source,val_time=value,t,k,candidate,why,perf_counter()-started
    globally_empty = globally_empty or lb > safe
    if not globally_empty:
        save(warm_plan(rt,inst["M"],inst["Q"]),"common_deterministic_dispatch")
        if initial is not None:
            save(initial, "user_supplied_replayed")
    raw=None
    if not globally_empty and ub!=lb:
        raw=[cert._windows(r,safe,inst["Q"]) for r in rt]
    num=0
    def cpstep(H,mode,cut,hint=None):
        nonlocal lb,num
        num+=1
        info,new,bound=cp_trial(inst,H,lb,folder/str(num),mode,cut,hint)
        if perf_counter()<=deadline:
            save(new,f"interval_{mode}_{num}")
            lb=max(lb,bound)
        history.append(dict(info,lower=lb,upper=ub,at=perf_counter()-started))
    if raw is not None:
        # Capacity implications use the absolute cap H on every route and every round trip.
        start=lb
        hi=ub if ub is not None else min(safe,lb+100)
        left,right=lb,hi
        while left<right and perf_counter()<deadline:
            H=(left+right)//2
            conflict=energy_check(inst,H,raw)
            if conflict is None:
                right=H
            else:
                left=H+1
                lb=max(lb,left)
                history.append(dict(method="absolute_window_energy",threshold=H,certificate=conflict,lower=lb))
        if ub!=lb and perf_counter()<deadline-1:
            cpstep(lb,"decision",min(deadline-.15,perf_counter()+2.5))
        if ub!=lb and lb<=safe and perf_counter()<deadline-.5:
            H=ub if ub is not None else min(safe,max(lb+20,ceil(lb*1.35)))
            cpstep(H,"optimize",min(deadline-.25,perf_counter()+min(15.,cap*.15)),plan)
        while ub!=lb and lb<=safe and perf_counter()<deadline-.7:
            H=ub-1 if ub is not None and ub-lb<=10 else (lb+ub)//2 if ub is not None else min(safe,max(lb+20,ceil(lb*1.6)))
            conflict=energy_check(inst,H,raw)
            if conflict:
                lb=H+1
                history.append(dict(method="absolute_window_energy",threshold=H,certificate=conflict,lower=lb))
                continue
            num+=1
            cutoff=min(deadline-.32,perf_counter()+min(12.,(deadline-perf_counter())*.25))
            try:
                info,new=compact_trial(inst,H,base,folder/str(num),cutoff)
            except TimeoutError:
                info,new=dict(method="compact",status="limit_build",threshold=H),None
            if perf_counter()<=deadline:
                save(new,f"compact_{num}")
                if info["status"]=="empty":
                    lb=max(lb,H+1)
            history.append(dict(info,lower=lb,upper=ub,at=perf_counter()-started))
            if ub==lb or perf_counter()>=deadline-.6: break
            if info["status"] not in ("empty","witness"):
                # Optimize over all improvements if an incumbent exists. Without one, the safe
                # global horizon covers an earliest realization of every feasible discrete plan.
                H=ub-1 if ub is not None else safe
                cpstep(H,"optimize",deadline-.15)
    globally_empty=globally_empty or lb>safe
    if ub is not None:
        assert not globally_empty and lb<=ub
    elapsed=perf_counter()-started
    return dict(status=result_status(lb,ub,globally_empty),
                lower=None if globally_empty else lb,upper=ub,plan=plan,truck_max=truck,last_event=last,plan_source=source,
                proof_seconds=elapsed if globally_empty or ub==lb else None,validation_seconds=val_time,
                elapsed_seconds=elapsed,time_limit=cap,history=history,global_safe_horizon=safe)
