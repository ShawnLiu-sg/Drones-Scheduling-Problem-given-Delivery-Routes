"""Run a scheduling algorithm on one published problem instance."""
import argparse
from contextlib import ExitStack
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

from scheduling import ablation, compact_flow, endurance, hybrid, liu_flow, makespan, relaxation
from scheduling.study import model, reference

METHODS = ("full","base","no_interval","no_capacity","liu","graph_a","price_d",
           "reference","fixed_policy","interval","compact","cbs","cbs_dp","pilot_milp","makespan",
           "endurance","release","price_original","price_filter","price_signature",
           "price_closure","price_all","price_core","price_adaptive","price_milp")

def input_problem(path):
    record=json.loads(path.read_bytes())
    problem=dict(record.get("instance",record))
    name=record.get("scenario",{}).get("id") or record.get("id") or problem.get("name",path.stem)
    horizon=problem.get("H")
    baseline=model.timetable(problem["routes"],[[] for _ in problem["routes"]])
    if baseline is not None:
        for route,events in zip(problem["routes"],baseline,strict=True):
            route["base"]=events["D"][-1]
    problem["H"]=makespan.global_horizon(problem["routes"]) if horizon is None else horizon
    return record,problem,name,baseline

def local_plan(path):
    if path is None:
        return None
    data=json.loads(path.read_bytes())
    if isinstance(data,list):
        return data
    plan=data.get("plan",data.get("routes"))
    if plan is None:
        raise ValueError("The plan file contains no physical witness.")
    return plan

def relaxation_bound(problem):
    return relaxation.lower_bound(problem["routes"])

def prescribed_start(problem):
    routes, policy = problem["routes"], problem["fixed_policy"]
    if len(policy) != len(routes):
        raise ValueError("The policy needs one parcel-to-stop list per truck.")
    for route, stops in zip(routes, policy):
        parcels = route["parcels"]
        if len(stops) != len(parcels) or any(
            type(stop) is not int or stop not in model.eligible(route, parcel)
            for stop, parcel in zip(stops, parcels)
        ):
            raise ValueError("The prescribed handoff list contains an ineligible or missing stop.")
    jobs = [(k, stop) for k, stops in enumerate(policy) for stop in sorted(set(stops))]
    machines = problem["M"]
    chains = [jobs[k::machines] for k in range(min(machines, len(jobs)))]
    return model.timetable(routes, policy, chains=chains)

def delay_solve(problem,method,cap,warm,start):
    lower=relaxation_bound(problem)
    if not isfinite(lower):
        return dict(status="infeasible",bound=None,objective=None,plan=None,origin="release_relaxation")
    if warm is None:
        warm=(prescribed_start(problem) if method=="fixed_policy" else
              relaxation.warm_plan(problem["routes"],problem["M"],problem["Q"]))
        if warm is not None:
            try: reference.validate(problem,warm)
            except ValueError: warm=None  # A legacy instance can impose a hard planning horizon.
    else:
        if method=="fixed_policy" and [p["assign"] for p in warm]!=problem["fixed_policy"]:
            raise ValueError("The supplied initial schedule uses different handoff sites.")
        reference.validate(problem,warm)
    upper=reference.validate(problem,warm) if warm is not None else None
    if warm is not None:
        tau=max(max(r["tau"]) for r in problem["routes"])
        # Any weakly improving solution has C_k <= base_k + total incumbent delay.
        # A replenishment return can be at most tau later than the truck's handoff.
        problem["H"]=min(problem["H"],max(r["base"] for r in problem["routes"])+upper+tau)
    problem["warm"]=warm
    if upper==lower:
        return dict(status="optimal",bound=lower,objective=upper,plan=warm,origin="common")
    cap=max(0,cap-(perf_counter()-start))
    if cap==0:
        return dict(status="limit",bound=lower,objective=upper,plan=warm)
    if method in ("full","base","no_interval","no_capacity"):
        # As in the timed experiment, check the cutoff during wide time-index builds.
        deadline=perf_counter()+cap
        original_clock=compact_flow.Clock
        def checked_clock(*args,**kwargs):
            if perf_counter()>=deadline: raise TimeoutError
            return original_clock(*args,**kwargs)
        def checked_ticks(seq):
            for i,tick in enumerate(seq):
                if i%2048==0 and perf_counter()>=deadline: raise TimeoutError
                yield tick
        def checked_range(*args):
            seq=range(*args)
            return checked_ticks(seq) if len(seq)>=16384 else seq
        def retain_initial(inst,_cap,*,threads):
            return dict(objective=upper,bound=0,routes=warm)
        compact_flow.range=checked_range
        try:
            with ExitStack() as patches:
                patches.enter_context(ablation.patch_attribute(compact_flow,"Clock",checked_clock))
                if method=="base":
                    patches.enter_context(ablation.patch_attribute(hybrid.interval_exact,"solve",retain_initial))
                result=ablation.solve(problem,cap,1,arm="no_capacity" if method=="base" else method)
        finally:
            del compact_flow.range
    elif method=="fixed_policy":
        result=reference.solve(problem,cap,1,fixed_assignments=problem["fixed_policy"])
        if result["routes"] is not None and [p["assign"] for p in result["routes"]]!=problem["fixed_policy"]:
            raise ArithmeticError("Returned handoffs differ from the prescribed policy.")
    elif method in ("liu","compact","interval","reference","graph_a","price_d","cbs","cbs_dp"):
        from scheduling import conflict_dp,conflict_search,interval_exact
        from scheduling.study import a_graph,d_price
        module={"liu":liu_flow,"compact":compact_flow,"interval":interval_exact,"reference":reference,
                "graph_a":a_graph,"price_d":d_price,"cbs":conflict_search,"cbs_dp":conflict_dp}[method]
        result=module.solve(problem,cap,1)
    elif method=="pilot_milp":
        from scheduling import pilot_milp
        result=pilot_milp.run(problem,cap)
    else:
        from scheduling import improved,engine
        result=(engine.run(problem,"capacity",cap) if method=="price_milp" else
                improved.solve(problem,method.removeprefix("price_"),cap))
    plan=result.pop("routes",None)
    if plan is not None and reference.validate(problem,plan)!=result["objective"]:
        raise ArithmeticError("Returned objective differs from the physical replay.")
    if result.get("bound") is not None:
        result["bound"]=max(result["bound"],lower)
    result["plan"]=plan
    return result

def endurance_solve(problem,record,name,initial,reference_total,cap,folder,start):
    routes=problem["routes"]
    target=sum(r["base"] for r in routes) if reference_total is None else reference_total
    state=dict(instance=problem,id=name,eta=record.get("eta",15),meta={},work_dir=folder,
               original_plan=initial,reference_total=target,plan_absolute=None,plan=None,
               physical="unknown",operational="unknown",absolute_lower=None,calls=[])
    state["meta"]["input_sha256"]=record["_input_sha256"]
    state["replay"]="missing"
    if initial is not None:
        # All non-endurance physics must already be valid for this fixed plan.
        unlimited=[dict(r,sorties=[dict(s,E=problem["H"]) for s in r["sorties"]]) for r in routes]
        reference.validate(dict(problem,routes=unlimited),initial)
        ok=endurance.save_plan(state,initial,"original_replay")
        state["replay"]="pass" if ok else "endurance_violation"
        state["violations"]=[dict(vehicle=k,sortie=j,away=u-l,limit=s["E"]) for k,(r,v) in enumerate(zip(routes,initial))
                             for j,(s,l,u) in enumerate(zip(r["sorties"],v["L"],v["U"])) if u-l>s["E"]]
    bound=relaxation_bound(problem)
    if not isfinite(bound):
        endurance.mark_infeasible(state,"release_time_relaxation")
    else:
        endurance.record_lower(state,sum(r["base"] for r in routes)+bound,"release_time_relaxation")
    state["precheck_seconds"]=perf_counter()-start
    endurance.resolve(state,cap)
    for field in ("instance","meta","work_dir","original_plan"):
        del state[field]
    return state

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method",choices=METHODS)
    parser.add_argument("input",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--seconds",type=float)
    parser.add_argument("--evidence-dir", type=Path, help="New local directory for makespan/endurance solver models and logs; defaults to temporary storage")
    parser.add_argument("--plan",type=Path,help="Locally saved plan from a compatible instance")
    parser.add_argument("--reference-total",type=int,help="Absolute truck-return sum for endurance feasibility")
    args=parser.parse_args()
    start=perf_counter()
    record,problem,name,bare=input_problem(args.input)
    record["_input_sha256"]=sha256(args.input.read_bytes()).hexdigest()
    cap=args.seconds if args.seconds is not None else (300 if record.get("scenario",{}).get("tier")=="large" else 120)
    if args.method=="fixed_policy" and "fixed_policy" not in problem:
        parser.error("fixed_policy requires an instance extracted with --fixed-policy.")
    if args.method!="fixed_policy" and "fixed_policy" in problem:
        parser.error("This instance prescribes handoffs. Use fixed_policy, or extract the unrestricted input.")
    plan=local_plan(args.plan)
    with TemporaryDirectory() as path:
        folder=Path(path)
        if args.evidence_dir is not None:
            folder=args.evidence_dir.resolve()
            folder.mkdir(parents=True,exist_ok=False)
        if bare is None or any(r["base"]>problem["H"] for r in problem["routes"]):
            result=dict(status="infeasible",plan=None,origin="bare_route_or_hard_horizon")
            if args.method=="endurance": result.update(physical="infeasible",operational="unattainable",replay="missing" if plan is None else "invalid")
        elif args.method=="makespan":
            result=makespan.solve(problem,max(0,cap-(perf_counter()-start)),folder,plan)
        elif args.method=="endurance":
            result=endurance_solve(problem,record,name,plan,args.reference_total,cap,folder,start)
        elif args.method=="release":
            from scheduling.release import evaluate
            if plan is None:
                parser.error("release requires --plan from a local run")
            reference.validate(problem,plan)
            result=evaluate(problem,plan,name)
        else:
            result=delay_solve(problem,args.method,cap,plan,start)
    result.update(input_id=record.get("id", name),arm=record.get("arm"),eta=record.get("eta"),method=args.method,scenario_id=name,input_sha256=record["_input_sha256"],elapsed_seconds=perf_counter()-start)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+"\n",encoding="utf-8")

if __name__=="__main__":
    main()
