"""Exact shared-fleet scheduling with A3 batch capacity and threshold search."""
from time import perf_counter

import gurobipy as gp

from scheduling import compact_flow as threshold
from scheduling import interval_exact
from scheduling import liu_flow as liu
from scheduling.a3_capacity import energy_cuts
from scheduling.certificate import _event_envelope


CAPACITY_ALGORITHM = "A3_4_dynamic"


def _individual_bounds(instance, **_kwargs):
    route_bounds = []
    for route in instance["routes"]:
        early, _, _, departure, _ = _event_envelope(route)
        route_bounds.append(max(0, early[departure[-1]] - route["base"]))
    bound = sum(route_bounds)
    return dict(bound=bound, route_bounds=route_bounds,
                other_bounds=[bound - delay for delay in route_bounds], intervals=[])


def solve(instance, time_limit=120, threads=1, *, progress=None, basic=False, trial=None):
    """Search delay budgets with optional A3 stages and a selected budget model."""
    start = perf_counter()
    routes = instance["routes"]
    n = len(routes)
    if trial is None:
        trial = threshold._trial
    if "H" not in instance or instance["M"] >= n:
        return threshold.solve(instance, time_limit, threads, progress=progress, trial=trial)

    deadline = start + time_limit
    best = instance.get("warm")
    upper = liu.validate(instance, best) if best is not None else None
    lower = 0
    stats = {"trials": []}

    def report():
        if progress is not None:
            progress(upper, lower)

    def finish(status):
        report()
        stats["seconds"] = perf_counter()-start
        return dict(status=status, objective=upper, bound=lower, routes=best,
                    origin="algorithm", stats=stats)

    if upper == 0:
        return finish("optimal")
    if best is None:
        return threshold.solve(instance, max(0, deadline - perf_counter()), threads,
                               progress=progress, trial=trial)

    if basic:
        energy = _individual_bounds(instance)
    else:
        fraction = 0.55 if 2 * instance["M"] <= n else 0.22
        interval_cap = min(max(0, deadline - perf_counter()), time_limit * fraction)
        primal = interval_exact.solve(instance, interval_cap, threads=threads)
        if primal["objective"] < upper:
            best, upper = primal["routes"], primal["objective"]
        lower = max(lower, primal["bound"])
        stats.update(interval_upper=upper, interval_lower=lower,
                     interval_seconds=perf_counter() - start)
        report()
        if upper == lower:
            return finish("optimal")
        slack = deadline - perf_counter()
        if slack <= 0:
            return finish("limit")
        # The caller's release-time bound guides the search for batches.
        # The capacity master certifies its own bound.
        gate_lower = max(lower, instance.get("_capacity_reference_lower", 0))
        energy = energy_cuts(instance, time_limit=min(1.8, slack * 0.25), threads=threads,
                             joint_bound=True, other_bounds=True, stage_bounds=(gate_lower, upper))
        if "stats" in energy:
            stats["energy_details"] = energy["stats"]
        stats.update(energy_lower=energy["bound"], energy_route_bounds=energy["route_bounds"],
                     energy_residual_bounds=energy["other_bounds"], energy_intervals=len(energy["intervals"]))
    lower = max(lower, energy["bound"])
    report()
    if lower == upper:
        return finish("optimal")
    assert lower<upper
    baselines=liu.schedules.timetable(routes,[[] for _ in routes])
    latest=[[min(liu._latest(r,p),len(r["tau"])-2) for p in r["parcels"]] for r in routes]
    setup=(baselines,latest,energy["route_bounds"])
    capacity = None if basic else energy
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag",0)
        env.start()
        while perf_counter() < deadline - 0.12:
            # Use all improving delays near closure and a smaller decision model for wide gaps.
            width = upper-lower
            target = upper-1 if width<=12 else lower+width//2
            try:
                state, witness, lb, info = trial(
                    env, instance, setup, target, deadline-0.10, threads, energy=capacity)
            except TimeoutError:
                break
            lower=max(lower,lb)
            if witness is not None:
                value=liu.validate(instance,witness)
                if value<upper:
                    best,upper=witness,value
            stats["trials"].append(dict(target=target,state=state,lower=lower,upper=upper,**info))
            report()
            assert lower<=upper
            if lower==upper:
                return finish("optimal")
            if state!="empty":
                break
    return finish("limit")
