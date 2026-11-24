"""Continuous event network for constructing the common starting schedules."""

def eligible(route, parcel):
    last = (parcel["dest"] if parcel["kind"] == "truck"
            else route["sorties"][parcel["dest"]]["l"])
    return [i for i, tau in enumerate(route["tau"]) if tau and i <= last]


def timetable(routes, assignments, *, fixed=None, chains=()):
    """Earliest solution with fixed parcel decisions and fixed starts or chains."""
    zero = ("time-zero",)
    ids = {zero: 0}
    edges = []

    def arc(a, b, w):
        edges.append((ids.setdefault(a, len(ids)), ids.setdefault(b, len(ids)), w))

    def eq(a, b, w):
        arc(a, b, w)
        arc(b, a, -w)

    for k, (rt, assign) in enumerate(zip(routes, assignments)):
        def event(kind, i):
            return (k, kind, i)
        n = len(rt["service"])
        eq(zero, event("A", 0), 0)
        eq(zero, event("D", 0), 0)
        for i in range(1, n):
            eq(event("D", i-1), event("A", i), rt["travel"][i-1])
            arc(event("A", i), event("D", i), rt["service"][i])
        for s, leg in enumerate(rt["sorties"]):
            for key, pos in (("L", leg["l"]), ("U", leg["b"])):
                arc(event("A", pos), event(key, s), 0)
                arc(event(key, s), event("D", pos), 0)
            arc(event("L", s), event("U", s), leg["F"])
            arc(event("U", s), event("L", s), -leg["E"])
            if s:
                arc(event("U", s-1), event("L", s), 0)
        for i in set(assign):
            tau = rt["tau"][i]
            arc(event("A", i), event("R", i), -tau)
            arc(event("R", i), event("D", i), tau)
            if fixed is not None:
                eq(zero, event("R", i), fixed[k][i])
        for p, i in zip(rt["parcels"], assign):
            tau = rt["tau"][i]
            arc(zero, event("R", i), p["r"])
            if p["kind"] == "sortie":
                arc(event("R", i), event("L", p["dest"]), tau)
            elif i == p["dest"]:
                arc(event("R", i), event("D", i), tau + rt["service"][i])
    for chain in chains:
        for (k, i), (l, j) in zip(chain, chain[1:]):
            arc((k, "R", i), (l, "R", j), 2*routes[k]["tau"][i])

    times = [0]*len(ids)
    for _ in times:
        changed = False
        for a, b, w in edges:
            if times[a] + w > times[b]:
                times[b] = times[a] + w
                changed = True
        if not changed:
            break
    else:
        return None
    out = []
    for k, (rt, assign) in enumerate(zip(routes, assignments)):
        sol = {}
        for key in ("A", "D", "L", "U"):
            size = len(rt["service"] if key in ("A", "D") else rt["sorties"])
            sol[key] = [times[ids[k,key,i]] for i in range(size)]
        sol["assign"] = list(assign)
        sol["starts"] = {i:times[ids[k,"R",i]] for i in sorted(set(assign))}
        out.append(sol)
    return out


def realize(route, assign, starts):
    result = timetable([route], [assign], fixed=[{int(i):t for i,t in starts.items()}])
    return result[0] if result else None


def baseline(route):
    rt = dict(route, parcels=[])
    sol = timetable([rt], [[]])
    if sol is None:
        raise ValueError("Delivery structure is infeasible without resupply")
    return sol[0]["D"][-1]
