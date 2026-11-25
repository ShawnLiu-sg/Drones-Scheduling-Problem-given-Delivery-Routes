"""Apply basic delay-budget search to Liu's replenishment network."""

from collections import defaultdict
from math import ceil
from time import perf_counter

import gurobipy as gp
from gurobipy import GRB

from scheduling import hybrid, liu_flow as liu


def _trial(env, instance, setup, budget, deadline, threads, *, energy=None):
    """Minimize total delay in the Liu network restricted to an integer budget."""
    if energy is not None:
        raise ValueError("The basic delay-budget search requires no capacity rows")

    baselines, latest, lower = setup
    total_lower = sum(lower)
    if total_lower > budget:
        return "empty", None, budget + 1, {}

    horizon = instance["H"]
    outgoing, incoming = defaultdict(list), defaultdict(list)
    stats = dict(handoff_states=0, departure_groups=0, parcel_variables=0,
                 depot_time_nodes=0, wait_arcs=0)
    build_start = perf_counter()
    with gp.Model("liu_delay_budget", env=env) as model:
        model.ModelSense = GRB.MINIMIZE
        model.Params.Threads = threads
        model.Params.Seed = 13
        model.Params.MIPGap = 0
        model.Params.MIPGapAbs = 0
        model.Params.FeasibilityTol = 1e-8
        model.Params.IntFeasTol = 1e-8
        model.Params.OptimalityTol = 1e-8
        teams = []
        for route, baseline, last, lb in zip(instance["routes"], baselines, latest, lower, strict=True):
            cap = min(budget - total_lower + lb, horizon - route["base"])
            if cap < lb:
                return "empty", None, budget + 1, {}
            team = liu._add_team(model, route, baseline, liu._road_tail(route), last,
                                 lb, cap, instance, outgoing, incoming, None, deadline, stats)
            if team is None:
                raise TimeoutError
            teams.append(team)

        if not liu._add_fleet_flow(model, outgoing, incoming, horizon, instance["M"], deadline, stats):
            raise TimeoutError
        model.update()
        # Liu's original objective sums the per-team terminal delays.
        model.addLConstr(model.getObjective() <= budget)
        model.update()
        stats["build_seconds"] = perf_counter() - build_start
        remaining = deadline - perf_counter()
        if remaining <= 0:
            raise TimeoutError
        model.Params.TimeLimit = remaining
        model.optimize()
        stats.update(variables=model.NumVars, binary=model.NumBinVars,
                     rows=model.NumConstrs, nodes=model.NodeCount, runtime=model.Runtime,
                     solver_status=model.Status)
        if model.Status == GRB.INFEASIBLE:
            return "empty", None, budget + 1, stats

        bound = total_lower
        if liu._finite_bound(model.ObjBound):
            bound = max(bound, min(budget + 1, ceil(model.ObjBound - liu.TOL)))
        witness = None
        if model.SolCount:
            witness = liu._witness(instance, model, teams)
            realized = liu.validate(instance, witness)
            if realized > budget or realized > model.ObjVal + liu.TOL:
                raise ArithmeticError("Physical replay exceeds the Liu budget or objective")
        state = "solved" if model.Status == GRB.OPTIMAL else "limit"
        return state, witness, bound, stats


def solve(instance, time_limit=120, threads=1, *, progress=None):
    """Return a verified schedule and lower bound using Liu's budgeted MILPs."""
    return hybrid.solve(instance, time_limit, threads, progress=progress, basic=True, trial=_trial)
