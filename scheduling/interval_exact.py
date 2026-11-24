"""Fixed-route FSTSP solved by optional depot-flight intervals and cumulative capacity.

The native IBM CP Optimizer 12.10 executable reads a CPO model and returns its
solution and global bound. No time-indexed start binaries or drone identities.
"""
from decimal import Decimal
from math import ceil, gcd
from pathlib import Path
from functools import reduce
import re
import subprocess
import sys
import tempfile
from time import perf_counter

from scheduling.study import model as clock
from scheduling.study.reference import validate

import os
CP = os.environ.get("CP_OPTIMIZER", "cpoptimizer")


def _array(values):
    return "[" + ", ".join(str(v) for v in values) + "]"


def _payload_units(instance):
    numbers = [Decimal(str(instance["Q"]))]
    numbers += [Decimal(str(p["w"])) for rt in instance["routes"] for p in rt["parcels"]]
    scale = 10 ** max(0, max(-n.as_tuple().exponent for n in numbers))
    units = [int(n * scale) for n in numbers]
    divisor = reduce(gcd, units)
    return {n: u // divisor for n, u in zip(numbers, units)}


def _build(instance, warm, upper):
    routes = instance["routes"]
    horizon = instance["H"]
    baselines = clock.timetable(routes, [[] for _ in routes])
    if baselines is None:
        return None
    units = _payload_units(instance)
    cap = units[Decimal(str(instance["Q"]))]
    lines, initial = [], []
    intervals, terminals, team_stops = [], [], []

    def integer(name, lo, hi, hint):
        lines.append(f"{name} = intVar({lo}..{hi});")
        if hint is not None:
            initial.append(f"{name} = ({int(hint)});")
        return name

    def add(expression):
        lines.append(expression + ";")

    for k, (rt, base) in enumerate(zip(routes, baselines)):
        seed = warm[k] if warm else None
        starts = {int(i): int(t) for i, t in seed["starts"].items()} if seed else {}
        n = len(rt["service"])
        road = [0] * n
        for i in reversed(range(n-1)):
            road[i] = road[i+1] + rt["travel"][i] + rt["service"][i+1]
        end = min(horizon, rt["base"] + upper) if upper is not None else horizon
        caps = [end - tail for tail in road]
        if any(v < e for v, e in zip(caps, base["D"])):
            return None

        def events(label, high):
            return [integer(f"{label}{k}_{i}", base[label][i], hi,
                            seed[label][i] if seed else None)
                    for i, hi in enumerate(high)]
        a = events("A", [c - s for c, s in zip(caps,rt["service"])])
        d = events("D", caps)
        l = events("L", [caps[f["l"]] for f in rt["sorties"]])
        u = events("U", [caps[f["b"]] for f in rt["sorties"]])
        add(f"{a[0]} == 0")
        add(f"{d[0]} == 0")
        for i in range(1,n):
            add(f"{a[i]} == {d[i-1]} + {rt['travel'][i-1]}")
            add(f"{d[i]} >= {a[i]} + {rt['service'][i]}")
        for j, f in enumerate(rt["sorties"]):
            for event, i in [(l[j], f["l"]), (u[j], f["b"])]:
                add(f"{a[i]} <= {event}")
                add(f"{event} <= {d[i]}")
            add(f"{u[j]} >= {l[j]} + {f['F']}")
            add(f"{u[j]} <= {l[j]} + {f['E']}")
            if j:
                add(f"{l[j]} >= {u[j-1]}")
        terminals.append(d[-1])
        choices = [clock.eligible(rt,p) for p in rt["parcels"]]
        candidates = sorted(set(i for row in choices for i in row))
        stops, loads, times, receipts, presences = [], [], [], [], []
        for i in candidates:
            tau = rt["tau"][i]
            first = max(0, base["A"][i]-tau,
                        min(p["r"] for p, row in zip(rt["parcels"], choices) if i in row))
            last = min(horizon - 2*tau, caps[i]-tau)
            if first > last:
                continue
            flight = f"j{k}_{i}"
            load = f"v{k}_{i}"
            add(f"{flight} = intervalVar(optional, start={first}..{last}, size={2*tau})")
            add(f"{load} = intVar(0..{cap})")
            pres = f"presenceOf({flight})"
            start = f"startOf({flight},0)"
            receipt = f"({start}+{tau})"
            add(f"{pres} == ({load}>0)")
            add(f"{pres} <= ({a[i]} <= {receipt})")
            add(f"{pres} <= ({receipt} <= {d[i]})")
            if seed:
                hint = "absent" if i not in starts else f"present, start={starts[i]}, end={starts[i]+2*tau}, size={2*tau}"
                initial.append(f"{flight} = ({hint});")
            stops.append(i)
            loads.append(load)
            times.append(start)
            receipts.append(receipt)
            presences.append(pres)
            intervals.append(f"pulse({flight},1)")
        team_stops.append(stops)
        assigns, weights = [], []
        for j, (p, opts) in enumerate(zip(rt["parcels"], choices)):
            domain = [index for index,i in enumerate(stops) if i in opts and
                      p["r"] <= min(horizon-2*rt["tau"][i],caps[i]-rt["tau"][i])]
            w = units[Decimal(str(p["w"]))]
            if not domain or w > cap:
                return None
            choice = f"p{k}_{j}"
            add(f"{choice} = intVar(" + ",".join(map(str,domain)) + ")")
            add(f"element({_array(times)}, {choice}) >= {p['r']}")
            received = f"element({_array(receipts)}, {choice})"
            if p["kind"] == "sortie":
                add(f"{l[p['dest']]} >= {received}")
            else:
                add(f"{d[p['dest']]} >= {received}+{rt['service'][p['dest']]}")
            if seed:
                initial.append(f"{choice} = ({stops.index(seed['assign'][j])});")
            assigns.append(choice)
            weights.append(w)
        if assigns:
            add(f"pack({_array(loads)}, {_array(assigns)}, {_array(weights)})")
            add(" + ".join(presences) + f" >= {(sum(weights)+cap-1)//cap}")
    if intervals:
        add(" + ".join(intervals) + f" <= {instance['M']}")
    cost = " + ".join(terminals) + f" - {sum(rt['base'] for rt in routes)}"
    if upper is not None:
        add(cost + f" <= {upper}")
    add(f"minimize({cost})")
    if initial:
        lines.extend(["startingPoint {", *initial, "}"])
    return "\n".join(lines), team_stops


def solve(instance, time_limit=120, threads=1, *, progress=None):
    begin = perf_counter()
    deadline = begin + time_limit
    best = instance.get("warm")
    upper = validate(instance,best) if best else None
    lower = 0
    stats = {}

    def finish(status):
        if progress:
            progress(upper, lower)
        stats["seconds"] = perf_counter()-begin
        return dict(status=status, objective=upper, bound=lower, routes=best, stats=stats,
                    origin="algorithm")

    if progress:
        progress(upper,lower)
    if upper == 0:
        return finish("optimal")
    if "H" not in instance:
        from scheduling.relaxation import lower_bound
        if lower_bound(instance["routes"]) == float("inf"):
            lower = None
            return finish("infeasible")
        raise ValueError("H is required for a feasible parcel relaxation")
    built = _build(instance, best, upper)
    if built is None:
        lower = None
        return finish("infeasible")
    cpo, stops = built
    stats.update(intervals=sum(map(len,stops)), build_seconds=perf_counter()-begin)
    with tempfile.TemporaryDirectory(prefix="fstsp-interval-") as folder:
        Path(folder,"schedule.cpo").write_text(cpo,encoding="utf8")
        # Include model loading, interpreter startup, schedule retrieval and
        # independent validation in the requested wall-clock budget.
        remaining = deadline - perf_counter() - 0.30
        if remaining <= 0:
            return finish("limit")
        cmd = [str(CP), "-c", "read schedule.cpo", f"set Workers {threads}",
               f"set TimeLimit {remaining:.4f}", "set OptimalityTolerance 0",
               "set RelativeOptimalityTolerance 0", "optimize", "display solution"]
        completed = subprocess.run(cmd, cwd=folder, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        log = completed.stdout.decode("utf8", errors="replace")
    status = re.search(r"Best objective\s*:\s*(-?\d+(?:\.\d+)?)",log)
    bound = re.search(r"Best bound\s*:\s*(-?\d+(?:\.\d+)?)",log)
    if bound:
        lower = max(0,ceil(float(bound.group(1))-1e-6))
    if status:
        assignments = {(int(k),int(p)): int(j)
                       for k,p,j in re.findall(r"^p(\d+)_(\d+)\[(\d+)\]",log,re.M)}
        dispatch = {(int(k),int(i)): int(t) for k,i,t in
                    re.findall(r"^j(\d+)_(\d+)\[1:\s*(\d+)\s+--",log,re.M)}
        fixed, chosen = [], []
        for k, (rt, sites) in enumerate(zip(instance["routes"],stops)):
            row = [sites[assignments[k,p]] for p in range(len(rt["parcels"]))]
            chosen.append(row)
            fixed.append({i:dispatch[k,i] for i in set(row)})
        candidate = clock.timetable(instance["routes"],chosen,fixed=fixed)
        objective = validate(instance,candidate)
        if objective > float(status.group(1))+1e-6:
            raise ArithmeticError("Replayed schedule exceeds the CP objective")
        if upper is None or objective < upper:
            best, upper = candidate, objective
    if "Optimal solution found" in log:
        assert upper is not None and lower == upper
        return finish("optimal")
    if upper is not None and lower == upper:
        return finish("optimal")
    if re.search(r"^ ! Search completed, model is infeasible|^Model is infeasible|^ ! Search completed, model has no solution",log,re.M):
        assert upper is None
        lower = None
        return finish("infeasible")
    if "Search terminated" not in log and "Search stopped" not in log:
        raise RuntimeError("CP Optimizer did not return an interpretable result:\n" + log[-2500:])
    return finish("limit")
