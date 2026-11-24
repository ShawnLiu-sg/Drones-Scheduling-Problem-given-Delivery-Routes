"""Exact route projection shared by all parcel-to-stop patterns of one team.

Static event paths are summarized once. For pattern vertices u,v, each walk
from u into the static network and back to v has best lag
max_{a,b} [out(u,a) + D_static(a,b) + in(b,v)]. Zero and the terminal
are retained with identity connectors. Closing this (B+2)-vertex graph also
accounts for walks using several dispatch vertices. A positive cycle here or
in the static network is equivalent to infeasibility in the full network.
"""
from time import perf_counter
from types import SimpleNamespace

from scheduling.fast_network import fast_distances
from scheduling.study.validation_d.temporal import Network


class ProjectedRoute:
    def __init__(self, rt, horizon, flight_horizon=None, deadline=None):
        self.route = rt
        self.horizon = horizon
        self.flight_horizon = flight_horizon
        self.static = Network(dict(rt, parcels=[]), [], 1, horizon, flight_horizon)
        graph = self.static
        # Every path leaving a dispatch enters D_i, the assigned sortie's L_s,
        # or zero. End and zero are the two static vertices retained in the core.
        # Their distances to every vertex suffice for all future patterns.
        source_set = {graph.end}
        source_set.update(graph.ids["D", i] for i, tau in enumerate(rt["tau"]) if tau)
        source_set.update(graph.ids["L", p["dest"]] for p in rt["parcels"] if p["kind"] == "sortie")
        graph.starts = sorted(source_set - {graph.zero})
        self.base_dist = fast_distances(graph, allow_infeasible=True, deadline=deadline)

    def pattern(self, row, fleet, deadline=None):
        if deadline is not None and perf_counter() >= deadline:
            raise TimeoutError("projected pattern")
        rt, base = self.route, self.static
        stops = sorted(set(row))
        durations = [2*rt["tau"][stop] for stop in stops]
        count = len(stops)+2
        net = SimpleNamespace(zero=0, end=1, starts=list(range(2,count)),
                              stops=stops, duration=durations, horizon=self.horizon)
        if self.base_dist is None:
            return net, None
        release = {stop:0 for stop in stops}
        outgoing = {stop:{base.ids["D",stop]:rt["tau"][stop]} for stop in stops}
        for p,stop in zip(rt["parcels"],row):
            release[stop] = max(release[stop], p["r"])
            tau = rt["tau"][stop]
            if p["kind"] == "sortie":
                outgoing[stop][base.ids["L",p["dest"]]] = tau
            elif p["dest"] == stop:
                outgoing[stop][base.ids["D",stop]] = tau + rt["service"][stop]
        if self.flight_horizon is not None:
            for stop, duration in zip(stops,durations):
                outgoing[stop][base.zero] = duration-self.flight_horizon

        exit_points = [{base.zero:0}, {base.end:0}] + [outgoing[stop] for stop in stops]
        arrivals = [base.ids["A",s] for s in stops]
        d = []
        for v, exits in enumerate(exit_points):
            sources = [(self.base_dist[a], lag) for a, lag in exits.items()]
            at_zero = max(static_row[base.zero]+lag for static_row,lag in sources)
            at_end = max(static_row[base.end]+lag for static_row,lag in sources)
            limits = [at_zero,at_end]
            for stop,a in zip(stops,arrivals):
                at_truck = max(static_row[a]+lag for static_row,lag in sources)
                limits.append(max(at_truck-rt["tau"][stop], at_zero+release[stop]))
            limits[v] = max(0, limits[v])
            d.append(limits)
        if fleet is not None:
            for pos in range(len(stops)-fleet):
                a, b = 2+pos, 2+pos+fleet
                d[a][b] = max(d[a][b], durations[pos])

        for k, middle in enumerate(d):
            if deadline is not None and perf_counter() >= deadline:
                raise TimeoutError("projected pattern closure")
            for src, dist in enumerate(d):
                through = dist[k]
                for v, lag in enumerate(middle):
                    candidate = through + lag
                    if candidate > dist[v]:
                        dist[v] = candidate
                if dist[src] > 0:
                    return net, None
        return net, d
