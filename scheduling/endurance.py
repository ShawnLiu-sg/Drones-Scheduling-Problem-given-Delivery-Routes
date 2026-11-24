"""Physical and reference-return feasibility for one endurance configuration."""
import json
import re
import subprocess
from hashlib import sha256
from math import isfinite
from pathlib import Path
from time import perf_counter
from scheduling import interval_exact as intervals
from scheduling.relaxation import _route_bound, warm_plan
from scheduling.study import model, reference

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def write(path, item):
    path.write_text(json.dumps(item, ensure_ascii=False), encoding="utf8")

def total(plan):
    return sum(team["D"][-1] for team in plan)


def plan_value(instance, plan):
    if plan is None:
        return None
    try:
        reference.validate(instance, plan)
    except ValueError:
        return None
    return total(plan)


def save_plan(state, plan, rule, source=None):
    value = plan_value(state["instance"], plan)
    if value is None:
        return False
    assert state["physical"] != "infeasible"
    assert state["absolute_lower"] is None or value >= state["absolute_lower"]
    if state["plan_absolute"] is None or value < state["plan_absolute"]:
        state.update(plan=plan, plan_absolute=value, plan_source=dict(rule=rule, source=source))
    if state["physical"] != "feasible":
        state.update(physical="feasible", physical_evidence=dict(rule=rule, source=source))
    if value <= state["reference_total"] and state["operational"] != "attain":
        assert state["operational"] != "unattainable"
        state.update(operational="attain", operational_evidence=dict(rule=rule, source=source))
    return True


def mark_infeasible(state, rule, source=None):
    assert state["physical"] != "feasible" and state["operational"] != "attain"
    evidence = dict(rule=rule, source=source)
    state.update(physical="infeasible", physical_evidence=evidence,
                 operational="unattainable", operational_evidence=evidence)


def record_lower(state, value, rule, source=None):
    if state["plan_absolute"] is not None:
        assert value <= state["plan_absolute"]
    if state["absolute_lower"] is None or value > state["absolute_lower"]:
        state.update(absolute_lower=value, lower_evidence=dict(rule=rule, source=source))
    if value > state["reference_total"] and state["operational"] != "unattainable":
        assert state["operational"] != "attain"
        state.update(operational="unattainable", operational_evidence=dict(rule=rule, source=source))


def monotone(states):
    """A plan must replay on the receiving configuration before it is reused."""
    for left in states:
        for right in states:
            if left is right:
                continue
            same = left["meta"]["physical_key"] == right["meta"]["physical_key"]
            if same:
                assert [[f["E"] for f in r["sorties"]] for r in left["instance"]["routes"]] == [
                    [f["E"] for f in r["sorties"]] for r in right["instance"]["routes"]]
            if left["eta"] < right["eta"] or same:
                better = left["plan_absolute"] is not None and (
                    right["plan_absolute"] is None or left["plan_absolute"] < right["plan_absolute"])
                if better:
                    assert save_plan(right, left["plan"], "identical_endurance_replay" if same
                                     else "tighter_endurance_replay", left["id"])
            if left["eta"] > right["eta"] or same:
                if left["physical"] == "infeasible" and right["physical"] != "infeasible":
                    mark_infeasible(right, "identical_endurance_certificate" if same
                                    else "looser_endurance_infeasibility", left["id"])
                if left["absolute_lower"] is not None and right["physical"] != "infeasible":
                    record_lower(right, left["absolute_lower"], "identical_endurance_bound" if same
                                 else "looser_endurance_bound", left["id"])


def retime(state):
    old = state["original_plan"]
    if old is None:
        return
    routes = state["instance"]["routes"]
    assignments = [row["assign"] for row in old]
    starts = [{int(i): t for i, t in row["starts"].items()} for row in old]
    plan = model.timetable(routes, assignments, fixed=starts)
    if save_plan(state, plan, "retime_fixed_batches_and_dispatches"):
        return
    # Recover anonymous aircraft paths from the old flight intervals.
    free = [0] * state["instance"]["M"]
    lanes = [[] for _ in free]
    flights = sorted((t, k, i) for k, row in enumerate(starts) for i, t in row.items())
    for t, k, i in flights:
        lane = min(range(len(free)), key=free.__getitem__)
        assert free[lane] <= t
        lanes[lane].append((k, i))
        free[lane] = t + 2 * routes[k]["tau"][i]
    save_plan(state, model.timetable(routes, assignments, chains=lanes),
              "retime_fixed_batches_and_aircraft_sequences")


def read_cp_plan(log, instance, stops):
    assigned = {(int(k), int(p)): int(i)
                for k, p, i in re.findall(r"^p(\d+)_(\d+)\[(\d+)\]", log, re.M)}
    flights = {(int(k), int(i)): int(t)
               for k, i, t in re.findall(r"^j(\d+)_(\d+)\[1:\s*(\d+)\s+--", log, re.M)}
    if len(assigned) != sum(len(r["parcels"]) for r in instance["routes"]):
        return None
    choices = [[sites[assigned[k, j]] for j in range(len(rt["parcels"]))]
               for k, (rt, sites) in enumerate(zip(instance["routes"], stops, strict=True))]
    if any((k, i) not in flights for k, row in enumerate(choices) for i in row):
        return None
    starts = [{i: flights[k, i] for i in set(row)} for k, row in enumerate(choices)]
    return model.timetable(instance["routes"], choices, fixed=starts)


def exact_decision(state, phase, deadline):
    begun = perf_counter()
    if deadline - begun < 0.6:
        return
    inst = state["instance"]
    budget = state["reference_total"] - sum(rt["base"] for rt in inst["routes"]) \
        if phase == "operational" else None
    path = state["work_dir"] / f"{phase}_{len(state['calls'])+1}"
    path.mkdir(parents=True, exist_ok=True)
    evidence = str(path.relative_to(state["work_dir"]))
    call = dict(phase=phase, cap_delay=budget, horizon=inst["H"], evidence=evidence,
                status="limit")
    state["calls"].append(call)
    built = intervals._build(inst, None, budget)
    plan = None
    elapsed = perf_counter() - begun
    if built is None:
        if perf_counter() <= deadline:
            call["status"] = "infeasible"
            write(path / "certificate.json", dict(reason="deterministic_empty_event_or_parcel_domain",
                                                 input_sha256=state["meta"]["input_sha256"],
                                                 cap_delay=budget, H=inst["H"]))
    else:
        cpo, stops = built
        lines = cpo.splitlines()
        assert sum(line.startswith("minimize(") for line in lines) == 1
        file = path / "schedule.cpo"
        file.write_text("\n".join(line for line in lines if not line.startswith("minimize(")) + "\n",
                        encoding="utf8")
        remaining = deadline - perf_counter()
        call["model_sha256"] = digest(file)
        if remaining > .4:
            cmd = [str(intervals.CP), "-c", "read schedule.cpo", "set Workers 1",
                   f"set TimeLimit {max(.05, remaining-.35):.4f}", "optimize", "display solution"]
            try:
                proc = subprocess.run(cmd, cwd=path, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      timeout=remaining)
                log = proc.stdout.decode("utf8", errors="replace")
            except subprocess.TimeoutExpired as error:
                log = (error.stdout or b"").decode("utf8", errors="replace")
            now = perf_counter()
            if now <= deadline:
                plan = read_cp_plan(log, inst, stops)
                value = plan_value(inst, plan)
                if perf_counter() <= deadline:
                    if value is not None:
                        assert budget is None or value <= state["reference_total"]
                        call["status"] = "feasible"
                    elif re.search(r"(?mi)^\s*(?:!\s*)?(?:Search completed, model (?:is infeasible|has no solution)|Model is infeasible)", log):
                        call["status"] = "infeasible"
                    elif not re.search(r"(?i)Search (?:stopped|terminated|completed)|Time limit", log):
                        raise ArithmeticError((state["id"], phase, log[-1500:]))
            (path / "solver.log").write_text(log, encoding="utf8")
        elapsed = perf_counter() - begun
    call["elapsed_seconds"] = elapsed
    write(path / "decision.json", call)
    if call["status"] == "feasible":
        assert save_plan(state, plan, "exact_cp_witness_" + phase, evidence)
    elif call["status"] == "infeasible":
        if phase == "physical":
            mark_infeasible(state, "exact_cp_global_infeasible", evidence)
        else:
            record_lower(state, state["reference_total"] + 1, "exact_cp_reference_infeasible", evidence)


def resolve(state, cap):
    begun = perf_counter()
    deadline = begun + cap - state["precheck_seconds"]
    if state["physical"] == "unknown" and state["original_plan"]:
        retime(state)
    if state["physical"] == "unknown":
        inst = state["instance"]
        save_plan(state, warm_plan(inst["routes"], inst["M"], inst["Q"]), "new_warm_dispatch")
    if state["operational"] == "unknown":
        cutoff = deadline if state["physical"] != "unknown" else begun + cap * .25
        exact_decision(state, "operational", min(cutoff, deadline))
    if state["physical"] == "unknown":
        exact_decision(state, "physical", deadline)
    if state["operational"] == "unknown":
        exact_decision(state, "operational", deadline)
    state["resolve_seconds"] = perf_counter() - begun
