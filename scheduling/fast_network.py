"""Exact rows of a longest-lag closure consumed by Scheme D's sparse timing solver.

The event network uses constraints x[v] >= x[u] + lag(u,v).  Bellman-Ford
from a virtual zero-weight super-source yields feasible potentials, or detects
a positive cycle anywhere in the graph.  The transformed lengths
h[v]-h[u]-lag(u,v) are nonnegative.  One Dijkstra search per needed source
then recovers all its original maximum lags exactly.  Missing source rows
are intentionally absent: the ordinary full-matrix MILP uses its own closure.
"""
from heapq import heappop, heappush
from time import perf_counter


def fast_distances(net, allow_infeasible=False, deadline=None):
    """Return {source: row} for the origin and every active handoff start."""
    count = len(net.names)
    potential = [0] * count
    for _ in range(count):
        if deadline is not None and perf_counter() >= deadline:
            raise TimeoutError("projected closure potentials")
        changed = False
        for u, v, lag in net.arcs:
            proposal = potential[u] + lag
            if proposal > potential[v]:
                potential[v] = proposal
                changed = True
        if not changed:
            break
    else:
        if allow_infeasible:
            return None
        raise ValueError("Infeasible base schedule in benchmark")

    outgoing = [[] for _ in range(count)]
    for u, v, lag in net.arcs:
        outgoing[u].append((v, potential[v] - potential[u] - lag))

    rows = {}
    for source in (net.zero, *net.starts):
        if deadline is not None and perf_counter() >= deadline:
            raise TimeoutError("projected closure search")
        shortest = [float("inf")] * count
        shortest[source] = 0
        pending = [(0, source)]
        while pending:
            length, u = heappop(pending)
            if length != shortest[u]:
                continue
            for v, cost in outgoing[u]:
                proposal = length + cost
                if proposal < shortest[v]:
                    shortest[v] = proposal
                    heappush(pending, (proposal, v))
        rows[source] = [h - potential[source] - d
                        for h, d in zip(potential, shortest)]
    if deadline is not None and perf_counter() >= deadline:
        raise TimeoutError("projected closure search")
    return rows
