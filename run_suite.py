"""Run a method serially on one published manifest; write newly computed outputs locally."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def plan_file(directory, case_id):
    return directory / (sha256(case_id.encode("utf-8")).hexdigest() + ".json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=["main", "cvrp", "endurance", "makespan"])
    parser.add_argument("method", help="Method accepted by run_solver.py, e.g. full or liu")
    parser.add_argument("--out", type=Path, help="Default: local_runs/DATASET/METHOD")
    parser.add_argument("--plans", type=Path, help="Reference-plan directory from run_suite.py main full")
    parser.add_argument("--seconds", type=float, help="Optional per-case compute budget in seconds")
    parser.add_argument("--ablation", action="store_true", help="Only cases with fewer supply UAVs than trucks")
    parser.add_argument("--tier", choices=["small","medium","large"], help="Select a size tier")
    parser.add_argument("--first", type=int, help="Process the first N matching manifest records")
    args = parser.parse_args()
    if args.method == "release" and not args.plans:
        parser.error("release requires --plans pointing to locally computed reference plans")
    directory = ROOT / "datasets" / args.dataset
    records = json.loads((directory/"manifest.json").read_text(encoding="utf-8-sig"))["records"]
    if args.tier:
        records = [r for r in records if r["tier"] == args.tier]
    if args.ablation:
        records = [r for r in records if r["M"] < r["K"]]
    if args.first is not None:
        records = records[:args.first]
    target = args.out or ROOT / "local_runs" / args.dataset / args.method
    target.mkdir(parents=True,exist_ok=True)
    processed = 0
    missing = []
    for record in records:
        case = directory / record["file"]
        if sha256(case.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Input hash mismatch: {case}")
        output = plan_file(target, record["id"])
        command = [sys.executable, str(ROOT/"run_solver.py"), args.method, str(case),
                   "--output", str(output)]
        if args.seconds is not None:
            command.extend(["--seconds",str(args.seconds)])
        if args.plans:
            ref = plan_file(args.plans, record["scenario_id"])
            available = ref.exists() and json.loads(ref.read_text(encoding="utf-8-sig")).get("plan") is not None
            if available:
                command.extend(["--plan",str(ref)])
            elif args.method == "release":
                missing.append(record["id"])
                continue
        subprocess.run(command,check=True)
        processed += 1
    print(json.dumps({"computed":processed,"missing_local_reference":missing},ensure_ascii=False))


if __name__ == "__main__":
    main()
