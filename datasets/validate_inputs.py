"""Audit the four formal input panels using only the Python standard library."""
from collections import Counter
from hashlib import sha256
import json
from math import ceil, dist
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCENARIO = set("id route source_instance instance family tier N K J_by_route P ratio rho_pct alpha alpha_tenths M".split())
ROUTE = set("route_id station_customer_ids sortie_customer_ids travel service tau sorties parcels".split())
PARCEL = set("kind dest w r customer demand_site release_factor_numerator".split())
SIZES = {"main":216, "cvrp":24, "endurance":864, "makespan":432}

def load(p):
    return json.loads(p.read_bytes())

def verify_input(j, original, dataset):
    additional = {"eta"} if dataset == "endurance" else {"arm"} if dataset == "makespan" else set()
    assert set(j) == {"id","scenario","instance","weight_unit_kg"} | additional
    sc, instance = j["scenario"], j["instance"]
    assert set(sc) == SCENARIO and set(instance) == {"name","Q","M","routes"}
    assert sc["M"] == instance["M"] and instance["Q"] == 6 and j["weight_unit_kg"] == .5
    assert len(instance["routes"]) == sc["K"] and sum(sc["J_by_route"]) == sc["P"]
    xy = {n["id"]:(n["x"],n["y"]) for n in original["nodes"]}
    depot=original["depot"]
    route_customers=[]
    for r, old, count in zip(instance["routes"], original["routes"], sc["J_by_route"], strict=True):
        assert set(r)==ROUTE and r["route_id"]==old["route_id"]
        seq, ids=r["station_customer_ids"],r["sortie_customer_ids"]
        assert seq[0]==seq[-1]==depot and 0 not in r["service"][1:-1]
        assert len(r["service"])==len(r["tau"])==len(seq)
        assert r["tau"]==[ceil(dist(xy[depot],xy[n])/2) for n in seq]
        assert r["travel"]==[ceil(sum(abs(x-y) for x,y in zip(xy[a],xy[b]))) for a,b in zip(seq,seq[1:])]
        route_customers.extend(seq[1:-1]+ids)
        assert set(seq[1:-1]+ids)==set(old["truck"][1:-1]+[s[1] for s in old["sorties"]])
        if j.get("arm")=="T":
            assert ids==[] and r["sorties"]==[]
            assert [x for x in seq if x in set(old["truck"])] == old["truck"]
        else:
            assert seq == old["truck"]
            for s, (a,b,c), customer in zip(r["sorties"],old["sorties"],ids,strict=True):
                assert set(s)=={"l","b","F","E"} and customer==b
                l,h=s["l"],s["b"]
                assert (seq[l],seq[h])==(a,c)
                assert s["F"]==ceil(dist(xy[a],xy[b])/2)+1+ceil(dist(xy[b],xy[c])/2)
                length=max(s["F"],sum(r["travel"][l:h])+sum(r["service"][l+1:h]))
                assert s["E"]==length+(length*j.get("eta",15)+99)//100
        assert len(r["parcels"])==count
        for p in r["parcels"]:
            assert set(p)==PARCEL and isinstance(p["r"],int) and p["r"]>=0 and p["w"] in (1,2,3)
            assert p["kind"] in ("truck","sortie")
            i=p["dest"]
            assert p["customer"]==(seq[i] if p["kind"]=="truck" else ids[i])
            assert p["demand_site"]==(i if p["kind"]=="truck" else r["sorties"][i]["l"])
            assert 7500<=p["release_factor_numerator"]<=12500
    assert len(set(route_customers))==len(route_customers)==sc["N"]
    assert set(route_customers)|{depot}==xy.keys()

def run():
    original={(r["family"],r["route"]):r for r in load(HERE/"source_routes.json")}
    assert len(original)==18
    frozen={}
    for dataset,n in SIZES.items():
        root=HERE/dataset
        m=load(root/"manifest.json")
        assert set(m)=={"dataset","records"} and m["dataset"]==dataset
        rows=m["records"]
        assert len(rows)==len({r["id"] for r in rows})==n
        assert {p.name for p in (root/"inputs").glob("*.json")}=={Path(r["file"]).name for r in rows}
        counts=Counter()
        for row in rows:
            file=root/row["file"]
            assert sha256(file.read_bytes()).hexdigest()==row["sha256"]
            j=load(file); sc=j["scenario"]
            assert row["id"]==j["id"] and row["scenario_id"]==sc["id"]
            for key in "route family tier N K M P rho_pct alpha".split(): assert row[key]==sc[key]
            verify_input(j, original[sc["family"],sc["route"]],dataset)
            parcels=[[(p["customer"],p["r"],p["w"],p["release_factor_numerator"]) for p in r["parcels"]] for r in j["instance"]["routes"]]
            if dataset=="main": frozen[sc["id"]]=parcels
            elif dataset in ("endurance","makespan"):
                assert frozen[sc["id"]]==parcels
                counts[sc["id"]]+=1
        if counts:
            assert counts.keys()==frozen.keys() and set(counts.values())=={4 if dataset=="endurance" else 2}
        print(f"{dataset}: {n} input manifests and physical records checked.")
    protocol=load(HERE/"release_uncertainty/protocol.json")
    assert protocol["kappa"]==[.05,.10,.20] and len(protocol["residual_distributions"])==4

if __name__=="__main__":
    run()
