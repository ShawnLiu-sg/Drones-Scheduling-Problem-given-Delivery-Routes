"""Safe domain reduction for a fixed-route pricing MILP or pattern enumerator.

Each surviving station may still be incompatible with other parcel decisions.
The returned lists preserve the canonical order of the original domains.
"""
from math import inf
from time import perf_counter

from scheduling.study.model import eligible
from scheduling.study.validation_d.temporal import Network


def _from_origin(origin, arcs, n, deadline):
    """Longest lower lags from origin; None denotes an inconsistent network."""
    dist = [-inf] * n
    dist[origin] = 0
    for _ in range(n):
        if deadline is not None and perf_counter() >= deadline:
            raise TimeoutError("pricing domain bounds")
        changed = False
        for u, v, lag in arcs:
            value = dist[u] + lag
            if value > dist[v]:
                dist[v] = value
                changed = True
        if not changed:
            return dist
    return None


def filtered_choices(rt, horizon, flight_horizon, deadline=None):
    """Keep all assignment options present in any feasible bounded route column."""
    parcels = rt["parcels"]
    baseline = Network(dict(rt, parcels=[]), (), 1, horizon)
    ids = baseline.ids
    origin = baseline.zero
    arcs = list(baseline.arcs)
    tau = rt["tau"]
    service = rt["service"]
    # A parcel received earlier on the route still needs truck travel to get
    # to its demand site. Under radial dominance, that travel offsets any
    # shorter flight to the earlier stop.
    dominance = all(d >= abs(a - b) for d, a, b in zip(rt["travel"], tau, tau[1:]))
    if dominance:
        for p in parcels:
            if p["kind"] == "truck":
                i = p["dest"]
                arcs.append((origin, ids["D",i], p["r"] + tau[i] + service[i]))
            else:
                j = p["dest"]
                i = rt["sorties"][j]["l"]
                arcs.append((origin, ids["L",j], p["r"] + tau[i]))
    earliest = _from_origin(origin, arcs, len(ids), deadline)
    if earliest is None:
        return [[] for _ in parcels]
    # T_0 = 0 and T_0 >= T_event + longest_lag(event,0).
    reverse = _from_origin(origin, [(v,u,w) for u,v,w in arcs], len(ids), deadline)
    latest = [-value for value in reverse]
    domains = []
    for p in parcels:
        kept = []
        for i in eligible(rt, p):
            # Earliest physically possible package receipt at this station.
            h = max(earliest[ids["A",i]], p["r"] + tau[i])
            if flight_horizon is not None and h + tau[i] > flight_horizon:
                continue
            extra = service[i] if p["kind"] == "truck" and p["dest"] == i else 0
            if h + extra > latest[ids["D",i]]:
                continue
            if p["kind"] == "sortie" and h > latest[ids["L",p["dest"]]]:
                continue
            kept.append(i)
        domains.append(kept)
    return domains
