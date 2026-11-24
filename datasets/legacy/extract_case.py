"""Extract one historical scheduling input without running a solver."""
import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("series")
    parser.add_argument("original_label")
    parser.add_argument("destination", type=Path)
    parser.add_argument("--fixed-policy",action="store_true",
                        help="add the prescribed handoff indices for a fixed-policy comparison")
    args = parser.parse_args()
    with (HERE / "index.csv").open(encoding="utf8", newline="") as fp:
        matches = [r for r in csv.DictReader(fp) if
                   (r["series"],r["original_label"]) == (args.series,args.original_label)]
    if len(matches)!=1:
        parser.error("The series and original label must match one row of index.csv.")
    entry = matches[0]
    where = entry["physical_id"]
    if entry["source"] == "formal":
        record = json.loads((HERE.parents[1] / where).read_bytes())
    else:
        with (HERE / "instances.jsonl").open(encoding="utf8") as fp:
            record = next(obj for line in fp if (obj := json.loads(line))["id"] == where)
    if args.fixed_policy:
        routes=record["instance"]["routes"]
        if entry["policy_rule"]=="demand_site":
            policy=[[p["dest"] if p["kind"]=="truck" else r["sorties"][p["dest"]]["l"]
                     for p in r["parcels"]] for r in routes]
        elif entry["policy_rule"]=="alternating_fixed_stops":
            policy=[]
            for r in routes:
                launch=r["sorties"][1]["l"]
                policy.append([1,launch,1,launch,4,5,6,7][:len(r["parcels"])])
        else:
            parser.error("This experimental input has no prescribed fixed-policy variant.")
        record["instance"]["fixed_policy"]=policy
    with args.destination.open("x", encoding="utf8") as fp:
        json.dump(record, fp, ensure_ascii=False, indent=2)
        fp.write("\n")

if __name__ == "__main__":
    main()
