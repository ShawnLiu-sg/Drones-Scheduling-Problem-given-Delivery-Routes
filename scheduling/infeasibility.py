from dataclasses import dataclass
Event = tuple[str, int]
ZERO: Event = ("zero", 0)


@dataclass(frozen=True)
class Edge:
    start: Event
    end: Event
    lag: int
    reason: str
    parcel: int | None = None


def time_edges(route):
    travel = route["travel"]
    service = route["service"]
    tau = route["tau"]
    if not all(road >= abs(right - left)
               for road, left, right in zip(travel, tau[:-1], tau[1:], strict=True)):
        raise ValueError("The release-time relaxation needs radial ground-time dominance.")
    edges = [
        Edge(ZERO, ("A", 0), 0, "depot_fixed_arrive"),
        Edge(("A", 0), ZERO, 0, "depot_fixed_arrive"),
        Edge(ZERO, ("D", 0), 0, "depot_fixed_depart"),
        Edge(("D", 0), ZERO, 0, "depot_fixed_depart"),
    ]
    for i, road in enumerate(travel, start=1):
        edges.extend((
            Edge(("D", i - 1), ("A", i), road, "truck_travel"),
            Edge(("A", i), ("D", i - 1), -road, "truck_travel"),
            Edge(("A", i), ("D", i), service[i], "truck_service"),
        ))
    for s, sortie in enumerate(route["sorties"]):
        launch, recover = sortie["l"], sortie["b"]
        edges.extend((
            Edge(("A", launch), ("L", s), 0, "launch_after_truck_arrival"),
            Edge(("L", s), ("D", launch), 0, "launch_before_truck_departure"),
            Edge(("A", recover), ("U", s), 0, "recovery_after_truck_arrival"),
            Edge(("U", s), ("D", recover), 0, "recovery_before_truck_departure"),
            Edge(("L", s), ("U", s), sortie["F"], "minimum_sortie_duration"),
            Edge(("U", s), ("L", s), -sortie["E"], "maximum_sortie_duration"),
        ))
        if s:
            edges.append(Edge(("U", s - 1), ("L", s), 0, "sortie_sequence"))
    for parcel in route["parcels"]:
        destination = parcel["dest"]
        if parcel["kind"] == "truck":
            lag = parcel["r"] + tau[destination] + service[destination]
            edges.append(Edge(ZERO, ("D", destination), lag, "truck_parcel_release",
                              parcel.get("customer")))
        else:
            launch = route["sorties"][destination]["l"]
            edges.append(Edge(ZERO, ("L", destination), parcel["r"] + tau[launch],
                              "sortie_parcel_release", parcel.get("customer")))
    return edges


def positive_cycle(edges):
    """Find and replay a strict contradiction using Bellman–Ford in max form."""
    vertices = {event for edge in edges for event in (edge.start, edge.end)}
    values = {vertex: 0 for vertex in vertices}
    predecessor = {}
    for _ in vertices:
        updated = None
        for edge in edges:
            lower = values[edge.start] + edge.lag
            if lower > values[edge.end]:
                values[edge.end] = lower
                predecessor[edge.end] = edge
                updated = edge.end
        if updated is None:
            return None
    on_cycle = updated
    for _ in vertices:
        on_cycle = predecessor[on_cycle].start
    cycle = []
    event = on_cycle
    while True:
        edge = predecessor[event]
        cycle.append(edge)
        event = edge.start
        if event == on_cycle:
            break
    cycle.reverse()
    assert all(a.end == b.start for a, b in zip(cycle, cycle[1:] + cycle[:1], strict=True))
    assert sum(edge.lag for edge in cycle) > 0
    return cycle


def depot_sortie_certificate(route):
    """A simpler release-versus-recovery-window certificate, when available."""
    candidates = []
    for s, sortie in enumerate(route["sorties"]):
        if sortie["l"] != 0:
            continue
        recover = sortie["b"]
        for parcel in route["parcels"]:
            i = parcel["dest"]
            if parcel["kind"] != "truck" or i >= recover:
                continue
            road = sum(route["travel"][i:recover])
            service = sum(route["service"][i:recover])
            earliest = parcel["r"] + route["tau"][i] + road + service
            if earliest > sortie["E"]:
                candidates.append(dict(
                    sortie_index_0=s, launch_stop_index_0=0, launch_fixed_tick=0,
                    recovery_stop_index_0=recover, customer=parcel.get("customer"),
                    truck_service_stop_index_0=i, release_tick=parcel["r"],
                    depot_to_customer_flight=route["tau"][i],
                    truck_travel_to_recovery=road, service_through_recovery_arrival=service,
                    recovery_lower=earliest, endurance_limit=sortie["E"],
                    violation=earliest - sortie["E"],
                ))
    return max(candidates, key=lambda row: row["violation"], default=None)
