"""Compute paired management comparisons from locally generated full-algorithm runs."""
import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path

from scheduling.comparisons import pair_status, triple_bounds

ROOT = Path(__file__).resolve().parent
TRANSITIONS = (
    ("late_fraction", "rho_pct", 20, 10),
    ("release_06_08", "alpha", .8, .6),
    ("release_08_10", "alpha", 1., .8),
    ("fleet_1_2", "M", 1, 2),
    ("fleet_2_3", "M", 2, 3),
)
FACTORS = ("family", "route", "rho_pct", "alpha", "M")


def local_evidence(row, runs):
    name = sha256(row["id"].encode("utf-8")).hexdigest() + ".json"
    path = runs / name
    if not path.exists():
        return {"status": "missing", "lower": None, "upper": None}
    r = json.loads(path.read_text(encoding="utf-8-sig"))
    if r.get("method") != "full" or r.get("input_sha256") != row["sha256"]:
        raise ValueError(f"Expected the matching full-algorithm run for {row['id']}")
    upper = r.get("objective") if r.get("plan") is not None else None
    lower = r.get("bound")
    state = r["status"].lower()
    if state == "optimal" and (upper is None or lower != upper):
        raise ValueError(f"Optimal record has no matching bounds and witness: {row['id']}")
    # Zero is a valid lower bound on total delay relative to the route baselines.
    return {"status": state, "lower": max(0, lower or 0), "upper": upper}


def write_table(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def group_by(rows, fixed):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in fixed)].append(row)
    return groups


def compare(rows, runs):
    evidence = {r["id"]: local_evidence(r, runs) for r in rows}
    pairs, triples = [], []
    for name, factor, av, bv in TRANSITIONS:
        fixed = [k for k in FACTORS if k != factor]
        for _, group in sorted(group_by(rows, fixed).items()):
            by_value = {r[factor]: r for r in group}
            a, b = by_value[av], by_value[bv]
            ea, eb = evidence[a["id"]], evidence[b["id"]]
            missing = [r["id"] for r, e in ((a, ea), (b, eb)) if e["status"] == "missing"]
            values = ((None,)*5 + ("未定", "缺少本地运行记录")) if missing else pair_status(ea, eb)
            lo, hi, mlo, mhi, exact, status, reason = values
            pairs.append(dict(comparison=name, family=a["family"], route=a["route"], tier=a["tier"], factor=factor,
                              A=a["id"], B=b["id"], A_status=ea["status"], B_status=eb["status"],
                              A_lower=ea["lower"], A_upper=ea["upper"], B_lower=eb["lower"], B_upper=eb["upper"],
                              raw_lower=lo, raw_upper=hi, monotone_lower=mlo, monotone_upper=mhi,
                              exact_difference=exact, classification=status, evidence=reason,
                              missing_ids=json.dumps(missing,ensure_ascii=False)))
    for key, group in sorted(group_by(rows, FACTORS[:-1]).items()):
        fleet = {r["M"]:r for r in group}
        ids = [fleet[m]["id"] for m in (1,2,3)]
        bound = [evidence[k] for k in ids]
        missing = [k for k, item in zip(ids,bound) if item["status"]=="missing"]
        values = ((None,)*4 + ("无法判断", "缺少本地运行记录")) if missing else triple_bounds(bound)
        rl, rh, ml, mh, status, reason = values
        triples.append(dict(zip(FACTORS[:-1], key), tier=fleet[1]["tier"], M1=ids[0], M2=ids[1], M3=ids[2],
                            raw_lower=rl, raw_upper=rh, monotone_lower=ml, monotone_upper=mh,
                            classification=status, evidence=reason, missing_ids=json.dumps(missing,ensure_ascii=False)))
    return pairs, triples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=ROOT/"local_runs/main/full")
    parser.add_argument("--out", type=Path, default=ROOT/"local_runs/management")
    args = parser.parse_args()
    records = json.loads((ROOT/"datasets/main/manifest.json").read_text(encoding="utf-8-sig"))["records"]
    pairs, triples = compare(records, args.runs)
    args.out.mkdir(exist_ok=True,parents=True)
    write_table(args.out/"single_factor_pairs.csv",pairs)
    write_table(args.out/"fleet_marginal_triples.csv",triples)
    counts = Counter((r["comparison"],r["classification"],bool(json.loads(r["missing_ids"]))) for r in pairs)
    counts.update(("fleet_marginal",r["classification"],bool(json.loads(r["missing_ids"]))) for r in triples)
    write_table(args.out/"counts.csv", [dict(comparison=c,classification=s,has_missing_input=missing,count=n)
                                        for (c,s,missing),n in sorted(counts.items())])
    print(f"Wrote {len(pairs)} comparisons and {len(triples)} triples to {args.out}")

if __name__ == "__main__":
    main()
