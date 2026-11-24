# Formal scheduling inputs

Each JSON file is one scheduling model input. The four manifests index 216 main scenarios (experiments 1, 2 and 7), 24 fixed CVRPLIB truck scenarios (experiment 3), 432 paired truck–drone and truck-only configurations (experiment 4) and 864 endurance configurations (experiment 5). The ablation panel uses the 120 main scenarios with M < K. Main, F-arm and 15% records intentionally reuse the same delivery model for matched comparisons. Experiment 6 samples independent parcel-release errors on main scenarios given a fixed schedule supplied at run time. Its seed, distributions, transformation and strengths are in `release_uncertainty/protocol.json`.

| Directory | Scenario key | Contents |
| --- | --- | --- |
| `main` | original scenario ID | Frozen routes, onboard delivery sorties, release times and shared fleet |
| `cvrp` | original scenario ID | Fixed all-truck routes with shared replenishment |
| `makespan` | original ID plus F or T | F uses the original truck and onboard-drone plan, T inserts drone-served customers into their original truck route |
| `endurance` | original ID plus e00, e05, e15 or e30 | Modified maximum away-from-truck duration for each onboard delivery sortie |

Each dataset has `manifest.json` with one record per input. `file` is relative to its manifest and `sha256` covers the JSON bytes. The top-level `scenario` describes the delivery input while `instance` holds the physical model. The input and reference clock use integer time ticks. Solvers compute their own route baselines, time horizon and any initial feasible schedule.

- `Q` and parcel `w` use half-kilogram units. The 3 kg shared-drone payload is encoded as `Q=6`. A parcel with `w=3` weighs 1.5 kg. `M` is the shared replenishment fleet size.
- `station_customer_ids` is the truck's ordered list of node IDs, including the repeated depot. Local station and sortie indices are zero-based. `travel` gives consecutive ground travel times, `service` gives each stop duration and `tau` gives the depot-to-station replenishment flight times. `sortie_customer_ids` identifies the customer assigned to each onboard sortie, in sortie order.
- Each onboard sortie has launch and recovery station indices `l,b`, minimum flight plus customer service time `F`, and maximum allowed elapsed time away from the truck `E`. Let T be the larger of F and the truck travel from l to b plus service at internal stops. For an endurance factor e expressed as an integer percent, E = T + (e T + 99) // 100. The main factor is 15%.
- A late parcel's absolute release time is `r`. Its `dest` is the station index for kind `truck`, or the sortie index for kind `sortie`. `demand_site` is the associated truck service or launch station. The parcel's frozen release jitter is `release_factor_numerator / 10000`. Paired scenarios retain the same customer identity, absolute release and weight when only fleet, endurance or delivery mode changes.
- `source_routes.json` contains the input node coordinates and fixed delivery structures for the 12 Ren et al. routes and six CVRPLIB routes. Ren data are identified by family and instance; see [Ren et al. (2024)](https://doi.org/10.1016/j.ejor.2024.05.031). The six CVRPLIB entries retain the original .vrp customer numbering, with depot 1; source links are recorded per instance.

Run `python datasets/validate_inputs.py` from the repository root with Python 3.10+. The command checks all file hashes, physical fields, customer coverage, route geometry, endurance rounding and consistency of the matched instances.
