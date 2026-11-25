"""Single-component variants of the Chapter 4 hybrid threshold algorithm.

These adapters retain the original hybrid solve's time budget, stopping
conditions and delay-budget selection. Each worker runs one solve serially.
"""
from contextlib import contextmanager
from time import perf_counter

from scheduling import hybrid


_full_solve = hybrid.solve


@contextmanager
def patch_attribute(owner, name, replacement):
    original = getattr(owner, name)
    setattr(owner, name, replacement)
    try:
        yield
    finally:
        setattr(owner, name, original)


def solve(instance, time_limit=120, threads=1, *, progress=None, arm):
    if arm not in {"full", "no_interval", "no_capacity"}:
        raise ValueError(f"Unknown solver variant: {arm}")
    if arm == "full" or instance["M"] >= len(instance["routes"]):
        return _full_solve(instance, time_limit, threads, progress=progress)

    begin = perf_counter()
    stages = {}

    if arm == "no_interval":
        original_energy = hybrid.energy_cuts

        def retain_incumbent(problem, _time_limit, *, threads):
            warm = problem["warm"]
            upper = sum(schedule["D"][-1] - route["base"]
                        for route, schedule in zip(problem["routes"], warm, strict=True))
            return dict(objective=upper, bound=0, routes=warm)

        def timed_energy(*args, **kwargs):
            start = perf_counter()
            result = original_energy(*args, **kwargs)
            finish = perf_counter()
            stages["capacity_stage"] = dict(start=start-begin, end=finish-begin,
                                            seconds=finish-start, lower=result["bound"],
                                            intervals=len(result["intervals"]),
                                            route_bounds=result["route_bounds"],
                                            other_bounds=result["other_bounds"])
            return result

        with patch_attribute(hybrid.interval_exact, "solve", retain_incumbent), \
             patch_attribute(hybrid, "energy_cuts", timed_energy):
            result = _full_solve(instance, time_limit, threads, progress=progress)
        omitted_prefix = "interval_"

    else:
        original_interval = hybrid.interval_exact.solve
        original_trial = hybrid.threshold._trial

        def timed_interval(*args, **kwargs):
            start = perf_counter()
            result = original_interval(*args, **kwargs)
            finish = perf_counter()
            stages["interval_stage"] = dict(start=start-begin, end=finish-begin,
                                            seconds=finish-start,
                                            lower=result["bound"], upper=result["objective"])
            return result

        def without_capacity(env, problem, setup, budget, deadline, threads, *, energy=None):
            return original_trial(env, problem, setup, budget, deadline, threads, energy=None)

        with patch_attribute(hybrid.interval_exact, "solve", timed_interval), \
             patch_attribute(hybrid, "energy_cuts", hybrid._individual_bounds), \
             patch_attribute(hybrid.threshold, "_trial", without_capacity):
            result = _full_solve(instance, time_limit, threads, progress=progress)
        omitted_prefix = "energy_"

    stats = result["stats"]
    for key in tuple(stats):
        if key.startswith(omitted_prefix):
            del stats[key]
    stats.update(stages)
    return result
