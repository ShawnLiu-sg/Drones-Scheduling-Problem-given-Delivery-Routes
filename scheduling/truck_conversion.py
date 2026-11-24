"""Deterministically insert delivery-drone customers into their assigned truck route."""
import math
def distance(xy, a, b):
    return math.ceil(sum(abs(x-y) for x, y in zip(xy[a], xy[b], strict=True)))

def insert_sorties(xy, original, sorties):
    route = list(original)
    remaining = [(index, 0, customer) for index, (_, customer, _) in enumerate(sorties)]
    trace = []
    while remaining:
        cost, sortie, position_in_sortie, chosen, edge = min(
            (distance(xy, a, j) + distance(xy, j, b) - distance(xy, a, b), s, t, j, i)
            for s,t,j in remaining for i,(a,b) in enumerate(zip(route, route[1:]))
        )
        trace.append(dict(step=len(trace)+1, sortie_index=sortie, sortie_customer_index=position_in_sortie,
                          customer=chosen, from_id=route[edge], to_id=route[edge+1], edge_index=edge, increment=cost))
        route.insert(edge+1, chosen)
        remaining.remove((sortie, position_in_sortie, chosen))
    return route, trace
