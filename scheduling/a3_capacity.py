"""Segmented interval-capacity certificates for the full threshold algorithm.

An option uses its segment's left endpoint as the delay cost and a certificate
at its right endpoint as its resource requirement. A3 refines the first
admissible segment when the layered relaxation stalls at the entry bound.
"""
from bisect import bisect_left, bisect_right
from collections import Counter
from time import perf_counter

if __package__:
    from . import certificate
    from .station_batch import batch_vector
else:
    import certificate
    from station_batch import batch_vector


def _intervals(instance, windows, max_extra):
    horizon = instance["H"]
    last = min(horizon, max((route["base"] + low + max_extra + max(route["tau"])
                             for route, (_, low) in zip(instance["routes"], windows)), default=0))
    intervals = [(a, a + width) for width in (1, 4, 8, 16, 32, 64)
                 for a in range(0, last - width + 1, max(1, width // 4))]
    candidates = set()
    for route, (raw, low) in zip(instance["routes"], windows):
        for parcel, row in zip(route["parcels"], raw):
            for site, (first, intercept, latest) in row.items():
                duration = 2 * route["tau"][site]
                end = min(latest, intercept + low + max_extra)
                if end < first:
                    continue
                for width in (duration, end - first + duration):
                    width = min(width, horizon)
                    for edge in (parcel["r"], first, end, first + duration, end + duration,
                                 first + duration // 2, end + duration // 2):
                        for start in (edge, edge - width):
                            a = max(0, min(start, horizon - width))
                            candidates.add((a, a + width))
    ordered = sorted(candidates.difference(intervals), key=lambda i: (i[1] - i[0], i[0]))
    count = min(96, len(ordered))
    intervals.extend(ordered[j * (len(ordered) - 1) // max(1, count - 1)] for j in range(count))
    return intervals


class _Grid:
    """Index the observation intervals by length for sparse overlap queries."""
    def __init__(self, intervals):
        by_width = {}
        for index, (a, b) in enumerate(intervals):
            by_width.setdefault(b - a, []).append((a, index))
        self.groups = []
        for width, pairs in by_width.items():
            starts, indices = zip(*sorted(pairs))
            self.groups.append((width, starts, indices))
        self.size = len(intervals)

    def positive(self, latest, finish, duration):
        for width, starts, indices in self.groups:
            left = bisect_right(starts, latest - width)
            right = bisect_left(starts, finish)
            for k in range(left, right):
                a = starts[k]
                yield indices[k], min(width, duration, a + width - latest, finish - a)


def _layered_vector(profile, grid, capacity):
    forced, pieces = profile
    fixed, free_volume, free_peak, weights = {}, {}, {}, {}
    for latest, finish, duration in forced:
        for j, part in grid.positive(latest, finish, duration):
            fixed[j] = fixed.get(j, 0) + part
    for latest, finish, duration, weight, free in pieces:
        for j, part in grid.positive(latest, finish, duration):
            histogram = weights.setdefault(j, {})
            histogram[part] = histogram.get(part, 0) + weight
            if free:
                free_volume[j] = free_volume.get(j, 0) + weight * part
                free_peak[j] = max(free_peak.get(j, 0), part)
    values = [0] * grid.size
    strict = max_gain = 0
    for j in fixed.keys() | weights.keys():
        histogram = weights.get(j, {})
        volume = sum(h * weight for h, weight in histogram.items())
        previous = max((volume + capacity - 1) // capacity, max(histogram, default=0),
                       fixed.get(j, 0) + max(free_peak.get(j, 0),
                                             (free_volume.get(j, 0) + capacity - 1) // capacity))
        load = sum(histogram.values())
        below = layered = 0
        for h in sorted(histogram):
            layered += (h - below) * ((load + capacity - 1) // capacity)
            load -= histogram[h]
            below = h
        values[j] = max(previous, layered)
        gain = values[j] - previous
        strict += int(gain > 0)
        max_gain = max(max_gain, gain)
    return tuple(values), strict, max_gain


def _build_tables(routes, windows, horizon, capacity, grid, max_extra, deadline, stats):
    neutral = (0,) * grid.size
    tables = [[(low, neutral)] for _, low in windows]
    first_segments = [None] * len(routes)
    completed = 0

    # Complete the fine prefix for each route before spending time on coarse
    # segments. A partial prefix retains the original neutral alternative.
    for k, (route, (raw, low)) in enumerate(zip(routes, windows)):
        if perf_counter() >= deadline:
            break
        fine_end = min(low + max_extra, horizon - route["base"])
        options = []
        previous = None
        first = None
        for d in range(low, fine_end + 1):
            if perf_counter() >= deadline:
                break
            profile = certificate._classify(route, raw, d, capacity)
            if profile is not None:
                if profile != previous:
                    vector, strict, gain = _layered_vector(profile, grid, capacity)
                stats["layer_strict_cells"] += strict
                stats["max_layer_gain"] = max(stats["max_layer_gain"], gain)
                if not options or vector != options[-1][1]:
                    options.append((d, vector))
                if first is None:
                    first = (d, d)
            previous = profile
        else:
            options.append((fine_end + 1, neutral))
            tables[k] = options
            first_segments[k] = first
            completed += 1
            continue
        break
    stats["neutral_fallback_routes"] = len(routes) - completed

    for k in range(completed):
        if perf_counter() >= deadline:
            break
        route = routes[k]
        raw, low = windows[k]
        high = horizon - route["base"]
        options = tables[k]
        for a in range(low + max_extra + 1, high + 1, 4):
            if perf_counter() >= deadline:
                break
            b = min(high, a + 3)
            profile = certificate._classify(route, raw, b, capacity)
            options.pop()  # Replace the old neutral suffix.
            if profile is not None:
                vector, strict, gain = _layered_vector(profile, grid, capacity)
                stats["layer_strict_cells"] += strict
                stats["max_layer_gain"] = max(stats["max_layer_gain"], gain)
                stats["nonzero_coarse_levels"] += int(any(vector))
                if not options or vector != options[-1][1]:
                    options.append((a, vector))
                if first_segments[k] is None:
                    first_segments[k] = (a, b)
            stats["coarse_levels_completed"] += 1
            options.append((b + 1, neutral))
        if len(options) > 1 and options[-1][0] == high + 1:
            options.pop()
    return tables, first_segments


def _master(tables, intervals, fleet, joint_bound, other_bounds, deadline, threads):
    local = [opts[0][0] for opts in tables]
    keep = [j for j, (a, b) in enumerate(intervals)
            if sum(opts[0][1][j] for opts in tables) > fleet * (b - a)]
    filtered = []
    for options in tables:
        levels = []
        for cost, vector in options:
            values = tuple(vector[j] for j in keep)
            if not levels or values != levels[-1][1]:
                levels.append((cost, values))
        filtered.append(levels)
    bound = sum(local)
    others = [bound - lb for lb in local]
    if keep and (joint_bound or other_bounds) and perf_counter() < deadline:
        capacities = [fleet * (intervals[j][1] - intervals[j][0]) for j in keep]
        with certificate.gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()
            if joint_bound:
                bound = certificate._master_bound(filtered, capacities, deadline, env, threads)
            if other_bounds:
                for k in range(len(tables)):
                    if perf_counter() >= deadline:
                        break
                    others[k] = certificate._master_bound(filtered[:k] + filtered[k + 1:],
                                                          capacities, deadline, env, threads)
    return local, keep, filtered, bound, others


def _refine_first_segments(tables, first_segments, routes, windows, capacity, intervals,
                           deadline, stats):
    outcomes = Counter()
    for k, (options, level) in enumerate(zip(tables, first_segments)):
        now = perf_counter()
        if now >= deadline:
            break
        if level is None:
            continue
        a, b = level
        floor = options[0][1]
        local_deadline = min(deadline, now + (deadline - now) / (len(routes) - k + 1))
        result = batch_vector(routes[k], windows[k][0], b, capacity, intervals,
                              deadline=local_deadline, certified_floor=floor)
        outcomes[result.status] += 1
        if result.status == "optimal":
            stronger = tuple(map(max, floor, result.values))
            if stronger == floor:
                continue
            gains = [new - old for old, new in zip(floor, stronger) if new > old]
            stats["batch_strict_cells"] += len(gains)
            stats["batch_max_gain"] = max(stats["batch_max_gain"], max(gains))
            if len(options) < 2 or options[1][0] > b + 1:
                options.insert(1, (b + 1, floor))
            options[0] = (a, stronger)
        elif result.status == "infeasible":
            stats["exact_batch_infeasible_levels"] += b - a + 1
            options.pop(0)
            if not options or options[0][0] > b + 1:
                options.insert(0, (b + 1, floor))
    stats["batch_status"] = dict(outcomes)


def energy_cuts(instance, time_limit=2.0, max_extra=12, threads=1, *,
                joint_bound=False, other_bounds=False, stage_bounds=(0, None)):
    """Compute A3/4/dynamic capacity rows, route, residual and joint bounds."""
    start = perf_counter()
    deadline = start + time_limit
    routes, horizon, capacity, fleet = (instance[key] for key in ("routes", "H", "Q", "M"))
    windows = [certificate._windows(route, horizon, capacity) for route in routes]
    intervals = _intervals(instance, windows, max_extra)
    independent = [low for _, low in windows]
    stats = dict(algorithm="A3_4_dynamic", independent_bounds=independent,
                 layer_strict_cells=0, max_layer_gain=0,
                 coarse_levels_completed=0, nonzero_coarse_levels=0, batch_status={},
                 batch_strict_cells=0, batch_max_gain=0, exact_batch_infeasible_levels=0)
    tables, first_segments = _build_tables(routes, windows, horizon, capacity, _Grid(intervals),
                                          max_extra, deadline, stats)
    local, keep, filtered, bound, others = _master(
        tables, intervals, fleet, joint_bound, other_bounds, deadline, threads)
    lower, upper = stage_bounds
    if upper is not None and max(lower, bound) >= upper:
        gate = "closed"
    elif bound > lower:
        gate = "A2_improved"
    elif perf_counter() >= deadline:
        gate = "stage_time_limit"
    else:
        gate = "A2_stalled"
    stats.update(batch_gate=gate, pre_batch_bound=bound, entry_lower=lower, entry_upper=upper)
    if gate == "A2_stalled":
        _refine_first_segments(tables, first_segments, routes, windows, capacity, intervals,
                               deadline, stats)
        if stats["batch_strict_cells"] or stats["exact_batch_infeasible_levels"]:
            a2_bound, a2_others = bound, others
            local, keep, filtered, bound, others = _master(
                tables, intervals, fleet, joint_bound, other_bounds, deadline, threads)
            bound = max(a2_bound, bound)
            others = list(map(max, others, a2_others))
    stats.update(excluded_below_route_bounds=sum(local) - sum(independent), energy_rows=len(keep),
                 choices=sum(map(len, filtered)), seconds=perf_counter() - start)
    return dict(intervals=[intervals[j] for j in keep], options=filtered, route_bounds=local,
                other_bounds=others, bound=bound, stats=stats)
