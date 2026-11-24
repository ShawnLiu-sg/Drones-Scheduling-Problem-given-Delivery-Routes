"""Matched, single-threaded solver arms for the formal Chapter 5 protocol."""
import math
import time
from collections import Counter

from scheduling.study import d_price, reference
from scheduling.study.validation_d.features import Strengthening
from scheduling.study.validation_d.temporal_global import EnumerationPricing


class Progress:
    """Monotone global certificates, timestamped from the start of model preparation."""

    def __init__(self, start):
        self.start = start
        self.trace = []
        self.upper = None
        self.lower = 0.0

    def __call__(self, upper, lower):
        changed = not self.trace
        if upper is not None and math.isfinite(upper) and abs(upper)<1e90:
            rounded = round(upper)
            if abs(upper-rounded)<1e-6:
                upper = rounded
            if self.upper is None or upper < self.upper:
                self.upper = upper
                changed = True
        if lower is not None and math.isfinite(lower) and abs(lower)<1e90 and lower > self.lower+1e-9:
            self.lower = lower
            changed = True
        if changed:
            self.trace.append(dict(seconds=time.perf_counter()-self.start,
                                   objective=self.upper,bound=self.lower))


class Routing:
    """Dispatch pricing only for a matching pricing domain and eligible node prices."""

    def __init__(self, I, hybrid, local):
        self.hybrid = hybrid
        self.capture = bool(I.get("_capture_duals"))
        self.enumerator = EnumerationPricing(I["Q"],I["M"] if local else None) if hybrid else None
        self.price_samples = []
        self.reasons = Counter()
        self.eligible_calls = 0
        self.eligible_team0 = 0
        self.milp_calls = 0

    def __call__(self,k,pricer,dual,prefix,phase,fixed,deadline,cut_costs):
        reasons = []
        if phase!=2:
            reasons.append("phase1")
        if fixed:
            reasons.append("branch_fixed")
        if cut_costs and any(cut_costs.values()):
            reasons.append("cut_price")
        if not reasons:
            self.eligible_calls += 1
            if k==0:
                self.eligible_team0 += 1
                if self.capture and self.eligible_team0 in (1,3,10):
                    self.price_samples.append(dict(ordinal=self.eligible_team0,team=k,phase=phase,
                        dual=dual,prefix=list(prefix),horizon=pricer.horizon,
                        flight_horizon=pricer.flight_horizon))
            if self.hybrid:
                return self.enumerator.price(pricer,dual,prefix,deadline)
        elif self.hybrid:
            self.reasons["+".join(reasons)] += 1
        self.milp_calls += 1
        return pricer.price(dual,prefix,phase,fixed,deadline,cut_costs)

    def stats(self):
        stats=dict(price_samples=self.price_samples,eligible_calls=self.eligible_calls,
                   fallback_reasons=dict(self.reasons),fallback_calls=sum(self.reasons.values()),
                   milp_pricing_calls=self.milp_calls)
        if self.enumerator is not None:
            stats.update({s:getattr(self.enumerator,s) for s in
                          ("enumerated","feasible","oracle_calls","complete_calls","fixed_solves",
                           "oracle_seconds","build_seconds")})
        return stats


def run(instance,arm,seconds):
    start=time.perf_counter()
    progress=Progress(start)
    if arm in ("direct","fixed_policy"):
        assignments=None
        if arm=="fixed_policy":
            assignments=instance.get("fixed_policy",instance.get("policy_assignments"))
            if assignments is None:
                assignments=[[route["sorties"][j]["l"] for j in (0,1,0,1)] for route in instance["routes"]]
        result=reference.solve(instance,max(0,seconds-(time.perf_counter()-start)),1,
                               fixed_assignments=assignments,progress=progress)
        result["stats"]["root_lp_exact"]=None
    elif arm in ("base","full","hybrid","local","capacity"):
        local=arm in ("full","local","capacity")
        strengthening=Strengthening(local=True,pairs=arm in ("full","capacity")) if local else None
        router=Routing(instance,arm in ("full","hybrid"),local)
        result=d_price.solve(instance,max(0,seconds-(time.perf_counter()-start)),1,
                             strengthening=strengthening,route_price=router,progress=progress)
        result["stats"].update(router.stats())
    else:
        raise ValueError("Unknown Chapter 5 arm: "+arm)
    elapsed=time.perf_counter()-start
    result["stats"]["wall_seconds"]=elapsed
    progress(result["objective"],result["bound"])
    if result["routes"] is not None:
        actual=reference.validate(instance,result["routes"])
        if actual!=result["objective"]:
            raise ArithmeticError(("Witness objective mismatch",actual,result["objective"]))
        if arm=="fixed_policy" and [r["assign"] for r in result["routes"]]!=assignments:
            raise ArithmeticError("Policy witness uses different handoff sites")
    if result["bound"] is not None and result["objective"] is not None and result["bound"]>result["objective"]+1e-6:
        raise ArithmeticError("Invalid global certificate")
    result["trace"]=progress.trace
    # Keep the independently replayed terminal record explicitly alongside all improvements.
    result["trace"].append(dict(seconds=elapsed,objective=result["objective"],bound=result["bound"]))
    return result
