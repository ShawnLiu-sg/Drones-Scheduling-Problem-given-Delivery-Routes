"""Release-time relaxation and feasible dispatch proposals for the public pilot."""
import math
import random
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

from scheduling.study import model, reference


def _site(rt, p):
    return p["dest"] if p["kind"] == "truck" else rt["sorties"][p["dest"]]["l"]


def _route_bound(rt):
    # Across any intervening truck segment, elapsed time dominates the
    # difference of the depot-to-site one-way flight times.
    if any(d < abs(a - b) for d, a, b in zip(rt["travel"], rt["tau"], rt["tau"][1:])):
        # Without radial dominance the no-resupply baseline is still a valid lower bound.
        return 0
    outgoing = defaultdict(list)
    root = ("zero", 0)
    def edge(a, b, duration):
        outgoing[a].append((b, duration))
        outgoing[b]
    def equality(a, b, duration):
        edge(a, b, duration)
        edge(b, a, -duration)
    equality(root, ("A",0), 0)
    equality(root, ("D",0), 0)
    for i in range(1, len(rt["service"])):
        equality(("D",i-1), ("A",i), rt["travel"][i-1])
        edge(("A",i), ("D",i), rt["service"][i])
    for s, f in enumerate(rt["sorties"]):
        for kind, pos in (("L",f["l"]), ("U",f["b"])):
            edge(("A",pos), (kind,s), 0)
            edge((kind,s), ("D",pos), 0)
        edge(("L",s), ("U",s), f["F"])
        edge(("U",s), ("L",s), -f["E"])
        if s:
            edge(("U",s-1), ("L",s), 0)
    for p in rt["parcels"]:
        i = _site(rt, p)
        if p["kind"] == "truck":
            edge(root, ("D",i), p["r"] + rt["tau"][i] + rt["service"][i])
        else:
            edge(root, ("L",p["dest"]), p["r"] + rt["tau"][i])

    # The super-source imposes all nonnegative event times. Keep its value
    # fixed to detect a route precluded by release and depot launch/endurance.
    times = {node:0 for node in outgoing}
    queue = deque(times)
    active = set(queue)
    depth = dict(times)
    while queue:
        a = queue.popleft()
        active.remove(a)
        for b, duration in outgoing[a]:
            earliest = times[a] + duration
            if earliest <= times[b]:
                continue
            if b == root or depth[a] + 1 >= len(times):
                return math.inf
            times[b] = earliest
            depth[b] = depth[a] + 1
            if b not in active:
                active.add(b)
                queue.append(b)
    return times["D", len(rt["service"])-1] - rt["base"]


def lower_bound(routes):
    """Fleet-relaxed lower bound on return delays. Infinity certifies infeasibility."""
    return sum(_route_bound(rt) for rt in routes)


def _lane_proposal(routes, rows, unconstrained, m, trial, rng):
    # Relative dispatches in each truck route are nondecreasing under radial
    # dominance. Trucks may change places between consecutive flights.
    pending = [deque(sorted(set(row))) for row in rows]
    lanes = [[] for _ in range(m)]
    free = [0] * m
    while any(pending):
        lane = min(range(m), key=free.__getitem__)
        options = []
        for k, sites in enumerate(pending):
            if sites:
                i = sites[0]
                arrival = unconstrained[k]["starts"][i]
                duration = 2 * routes[k]["tau"][i]
                if trial == 0:
                    score = arrival
                elif trial == 1:
                    score = max(free[lane], arrival) + duration
                else:
                    score = arrival + rng.uniform(-1,1) * duration
                options.append((score, k))
        _, k = min(options)
        i = pending[k].popleft()
        lanes[lane].append((k,i))
        free[lane] = max(free[lane], unconstrained[k]["starts"][i]) + 2*routes[k]["tau"][i]
    return lanes


def warm_plan(routes, m, Q=3):
    """Build and replay candidate schedules; failure returns None."""
    bound = lower_bound(routes)
    if not math.isfinite(bound):
        return None
    rng = random.Random(429719)
    demand = [[_site(rt,p) for p in rt["parcels"]] for rt in routes]
    best = None
    best_delay = math.inf
    for attempt in range(12):
        rows = [list(row) for row in demand]
        if attempt:
            for k, (row, rt) in enumerate(zip(rows, routes)):
                loads = Counter()
                for i, p in zip(row,rt["parcels"]):
                    loads[i] += p["w"]
                for j in rng.sample(range(len(row)), len(row)):
                    old, p = row[j], rt["parcels"][j]
                    earlier = [i for i in model.eligible(rt,p)
                               if i != old and loads[i] + p["w"] <= Q]
                    occupied = [i for i in earlier if loads[i]]
                    if occupied and rng.random() < .45:
                        new = max(occupied)
                    elif earlier and rng.random() < .18:
                        new = rng.choice(earlier)
                    else:
                        continue
                    loads[new] += p["w"]
                    loads[old] -= p["w"]
                    row[j] = new
        if any(max((sum(p["w"] for p,j in zip(rt["parcels"],row) if i==j)
                    for i in set(row)),default=0)>Q for row,rt in zip(rows,routes)):
            continue
        relaxed = model.timetable(routes, rows)
        if relaxed is None:
            continue
        for trial in range(4 if not attempt else 2):
            lanes = _lane_proposal(routes, rows, relaxed, m, trial, rng)
            plan = model.timetable(routes, rows, chains=lanes)
            if plan is None:
                continue
            delay = sum(w["D"][-1] - rt["base"] for rt,w in zip(routes,plan))
            if delay < best_delay:
                h = max(max(rt["base"] for rt in routes)+delay+max(max(rt["tau"]) for rt in routes),
                        max(max(w["D"]) for w in plan))
                I = dict(routes=routes, H=h, Q=Q, M=m)
                if reference.validate(I,plan) != delay:
                    raise ArithmeticError("Warm schedule objective failed replay")
                best, best_delay = plan, delay
                if delay == bound:
                    return best
    return best
