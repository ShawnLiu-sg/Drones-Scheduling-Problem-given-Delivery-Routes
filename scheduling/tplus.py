"""Build T+ truck tours before scheduling, using only route geometry and assignments.

Run from the repository root, for example:
    python -m scheduling.tplus rat99-m5 --output-dir routes_tplus
"""

import argparse
import json
from itertools import combinations, pairwise
from math import ceil
from pathlib import Path

import gurobipy as gp
from gurobipy import GRB

from .truck_conversion import distance, insert_sorties


def _shortest_cycle(edges):
    """Find a shortest component in a degree-two solution with symmetric edge keys."""
    adjacent = {}
    for first, second in edges:
        adjacent.setdefault(first, []).append(second)

    # Match the experiment's traversal order when several minimum-size cuts tie.
    remaining = set(iter(adjacent))
    shortest = None
    while remaining:
        node = list(remaining)[-1]
        cycle = []
        while node is not None:
            cycle.append(node)
            remaining.remove(node)
            node = next((v for v in reversed(adjacent[node]) if v in remaining), None)
        if shortest is None or len(cycle) < len(shortest):
            shortest = cycle
    return shortest


def _orient(cycle, original):
    depot = original[0]
    pivot = cycle.index(depot)
    forward = cycle[pivot + 1:] + cycle[:pivot]
    rank = {node: index for index, node in enumerate(original[1:-1])}

    def score(order):
        inversions = sum(rank[a] > rank[b] for a, b in combinations(order, 2))
        return inversions, order

    return [depot, *min(forward, forward[::-1], key=score), depot]


def optimize_tour(env, xy, original):
    """Minimize one truck's ground travel, retaining the original T tour on a tie."""
    original = list(original)
    nodes = sorted(original[:-1])
    if original[0] != original[-1] or len(nodes) != len(set(nodes)):
        raise ValueError("Expected a closed truck tour with distinct interior stops")
    initial = sum(distance(xy, a, b) for a, b in pairwise(original))
    sequence, travel, lower, status = original, initial, initial, "optimal"

    if len(nodes) > 2:
        costs = {edge: distance(xy, *edge) for edge in combinations(nodes, 2)}
        warm_edges = {tuple(sorted(edge)) for edge in pairwise(original)}
        with gp.Model("Tplus", env=env) as model:
            model.Params.Threads = 1
            model.Params.Seed = 13
            model.Params.LazyConstraints = 1
            model.Params.MIPGap = 0
            model.Params.MIPGapAbs = 0
            model.Params.WorkLimit = 35
            edges = model.addVars(costs.keys(), obj=costs, vtype=GRB.BINARY, name="edge")
            for edge, var in edges.items():
                var.Start = int(edge in warm_edges)
            edges.update({(b, a): var for (a, b), var in edges.items()})
            for a in nodes:
                model.addConstr(gp.quicksum(edges[a, b] for b in nodes if b != a) == 2)

            def cut_subtour(model, where):
                if where == GRB.Callback.MIPSOL:
                    values = model.cbGetSolution(edges)
                    cycle = _shortest_cycle(edge for edge, value in values.items() if value > 0.5)
                    if len(cycle) < len(nodes):
                        model.cbLazy(gp.quicksum(edges[e] for e in combinations(cycle, 2)) <= len(cycle) - 1)

            model.optimize(cut_subtour)
            if model.Status not in (GRB.OPTIMAL, GRB.WORK_LIMIT):
                raise RuntimeError(f"Unexpected Gurobi status {model.Status}")
            status = "optimal" if model.Status == GRB.OPTIMAL else "work_limit"
            lower = ceil(max(0, model.ObjBound) - 1e-6)
            if model.SolCount and round(model.ObjVal) < initial:
                cycle = _shortest_cycle(edge for edge, var in edges.items() if var.X > 0.5)
                if len(cycle) != len(nodes):
                    raise RuntimeError("T+ incumbent does not visit every assigned stop")
                sequence = _orient(cycle, original)
                travel = sum(distance(xy, a, b) for a, b in pairwise(sequence))

    return dict(original=original, sequence=sequence, original_travel=initial,
                plus_travel=travel, improvement=initial - travel, travel_lower=lower,
                status=status, exact=lower == travel)


def optimize_routes(env, source):
    """Build original T by minimum-increment insertion, then optimize each truck."""
    xy = {node["id"]: (node["x"], node["y"]) for node in source["nodes"]}
    routes = []
    for truck in source["routes"]:
        original, _ = insert_sorties(xy, truck["truck"], truck["sorties"])
        routes.append(dict(route_id=truck["route_id"], **optimize_tour(env, xy, original)))
    return dict(family=source["family"], route=source["route"], depot=source["depot"], routes=routes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("routes", nargs="+", help="source route names, e.g. rat99-m5")
    default_catalog = Path(__file__).resolve().parents[1] / "datasets" / "source_routes.json"
    parser.add_argument("--catalog", type=Path, default=default_catalog)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sources = {s["route"]: s for s in json.loads(args.catalog.read_text(encoding="utf-8"))}
    unknown = set(args.routes) - sources.keys()
    if unknown:
        parser.error(f"Unknown source routes: {', '.join(sorted(unknown))}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with gp.Env(empty=True) as env:
        env.setParam("LogToConsole", 0)
        env.start()
        for name in dict.fromkeys(args.routes):
            result = optimize_routes(env, sources[name])
            path = args.output_dir / f"{name}.json"
            path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
