"""Exact, bounded implicit station-batch search for one route/delay cap.

The result is a coordinate-wise minimum. Different intervals may attain their
minima on different assignments. Interrupted trees export no objective value.
"""
from dataclasses import dataclass
from time import perf_counter


@dataclass(frozen=True)
class BatchResult:
    status: str
    values: tuple | None
    nodes: int
    assignments: int


def batch_vector(route, raw, delay, capacity, intervals, *, deadline, max_nodes=4096,
                 certified_floor=None):
    """Use raw _windows rows and an optional proved lower bound for this local model."""
    if perf_counter() >= deadline:
        return BatchResult("time_limit", None, 0, 0)
    items = []
    for parcel, row in zip(route["parcels"], raw, strict=True):
        choices = tuple(sorted((i, f, min(t + delay, cap))
                               for i, (f, t, cap) in row.items()
                               if f <= min(t + delay, cap)))
        if not choices or parcel["w"] > capacity:
            return BatchResult("infeasible", None, 0, 0)
        items.append((parcel["w"], choices))
    items.sort(key=lambda x: (len(x[1]), -x[0], x[1]))
    active = {}
    best = None
    floor = tuple(certified_floor) if certified_floor is not None else (0,) * len(intervals)
    nodes = leaves = 0
    stop = None

    def search(j, previous_site):
        nonlocal best, nodes, leaves, stop
        if perf_counter() >= deadline:
            stop = "time_limit"
            return
        if nodes >= max_nodes:
            stop = "node_limit"
            return
        nodes += 1
        if j == len(items):
            leaves += 1
            values = [0] * len(intervals)
            for i, (_, left, right) in active.items():
                duration = 2 * route["tau"][i]
                for k, (a, b) in enumerate(intervals):
                    if k % 128 == 0 and perf_counter() >= deadline:
                        stop = "time_limit"
                        return
                    values[k] += max(0, min(b-a, duration, b-right, left+duration-a))
            best = values if best is None else list(map(min, best, values))
            if tuple(best) == floor:
                stop = "optimal"
            return
        weight, options = items[j]
        interchangeable = j > 0 and items[j] == items[j-1]
        for site, left, right in options:
            # Permuting packages with identical weight and every window
            # preserves the occupied stations and all interval objectives.
            if interchangeable and site < previous_site:
                continue
            original = active.get(site)
            load = weight
            if original is not None:
                used, first, last = original
                load += used
                left, right = max(left, first), min(right, last)
            if load > capacity or left > right:
                continue
            active[site] = (load, left, right)
            search(j+1, site)
            if stop is not None:
                return
            if original is None:
                del active[site]
            else:
                active[site] = original
    search(0, None)
    if stop is not None and stop != "optimal":
        return BatchResult(stop, None, nodes, leaves)
    if best is None:
        return BatchResult("infeasible", None, nodes, leaves)
    return BatchResult("optimal", tuple(best), nodes, leaves)
