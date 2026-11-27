# Exact scheduling of shared resupply drones for truck–drone tandems on fixed routes

This repository accompanies the paper of the same title. It contains scheduling algorithms and benchmark instances for delivery teams that follow fixed truck routes and delivery-drone sorties. A shared fleet of depot-based **resupply drones** carries **delayed-available (DA) parcels** to the teams after dispatch. The scheduler chooses parcel batches, rendezvous stops, and resupply flight times while enforcing parcel availability, payload and fleet capacity, truck–drone synchronization, and delivery-drone endurance.

The primary objective is to minimize **total return delay**, the sum of the delivery teams' return delays relative to their no-resupply baselines. The proposed exact time-threshold algorithm combines interval search, a layered load bound, and a **same-stop batch bound (SBB)**. The two main benchmarks use the time-expanded resupply network of Liu et al. (2024). **Liu** solves it directly, and **Liu + threshold** applies time-threshold budget search to the same model.

## Paper experiments

| Section | Experiment | Design and main finding |
| --- | --- | --- |
| §5.2 | Algorithm comparison on fixed truck–drone routes | Compare the proposed algorithm, Liu, and Liu + threshold on 216 scenarios. They prove optimality in 180, 172, and 173 scenarios, respectively. The proposed algorithm attains a 7.26% mean relative optimality gap and reduces mean capped time to proof by 24.3% versus Liu and 23.0% versus Liu + threshold. |
| §5.3 | Algorithm ablation | Compare four configurations on 120 scenarios with fewer resupply drones than delivery teams. Interval search accounts for most of the additional proofs. Adding SBB-based capacity strengthening to interval search reduces the mean gap from 14.65% to 10.54%. |
| §5.4 | Algorithm comparison on fixed truck-only routes | Compare the proposed algorithm and Liu on 24 scenarios from six CVRPLIB route sets. The methods prove 16 and 14 optima and attain mean gaps of 17.52% and 33.31%, respectively. |
| §5.5 | Makespan with matched customer assignments | Compare truck–drone tandems (F), insertion-based truck-only routes (T), and truck-only routes optimized within each original customer assignment (T+). With sufficient resupply capacity, F is faster than both truck-only configurations in all 72 matched settings. In the eight one-drone reversals, a second drone restores F's advantage in all eight against T and seven against T+. The remaining reversal against T+ requires three drones. |
| §5.6 | Delivery-drone endurance margins | Evaluate 216 scenarios at 0%, 5%, 15%, and 30% endurance margins. Starting from the 15% reference schedules, rescheduling restores feasibility for 132 of 140 invalidated schedules at 5% and 164 of 172 at 0%. |
| §5.7 | DA parcel availability-time perturbations | Use 208 verified reference schedules, four error distributions, and perturbation intensities of 5%, 10%, and 20%. With batches, rendezvous stops, and flight sequences fixed and waiting allowed, feasibility exceeds 96% at 5% intensity and ranges from 92.62% to 93.74% at 20%. |
| §5.8 | Sensitivity to resupply capacity and DA parcel parameters | Vary the DA parcel share, availability-time coefficient, or resupply fleet size one factor at a time. Adding a second and a third resupply drone reduces total return delay in 54 and 27 of 72 matched comparisons, respectively. |

In §5.2, mean relative optimality gaps use the 206 common scenarios where all three methods obtain finite lower and upper bounds. Mean capped time to proof uses the 156 independently solved scenarios after 60 common preprocessing proofs. In §5.5, makespan ends when both the last delivery truck and the last active resupply drone have returned to the depot.

## Run locally

Run the commands from the repository root. The experiments used Python, Gurobi 12.0.1, and IBM ILOG CP Optimizer 12.10. Use Python 3.11 and a working Gurobi license. For methods that invoke CP Optimizer, add `cpoptimizer` to `PATH` or set `CP_OPTIMIZER` to its executable path.

```sh
python -m pip install -r requirements.txt
python datasets/validate_inputs.py
python run_suite.py main full --tier small --first 1
python run_suite.py main liu --tier small --first 1
python run_suite.py main liu_budget --tier small --first 1
```

Omit `--first` and `--tier` to run the full manifest. The default per-case time limits are 120 seconds for small and medium cases and 300 seconds for large cases. Results go to `local_runs/<dataset>/<method>/` by default. Use `python run_solver.py --help` for single-instance options.

| Section | Repository entry points |
| --- | --- |
| §5.2 | Run `python run_suite.py main <method>` with `full` for the proposed algorithm, `liu` for Liu, and `liu_budget` for Liu + threshold. |
| §5.3 | Run `python run_suite.py main <method> --ablation --out local_runs/ablation/<method>` with `base` for the basic time-threshold method, `no_capacity` for interval search only, `no_interval` for SBB-based capacity strengthening only, and `full` for both components. |
| §5.4 | Run `python run_suite.py cvrp full` and `python run_suite.py cvrp liu`. |
| §5.5 | Run `python run_suite.py makespan makespan` on the 216 F/T pairs in the published makespan manifest. Generate T+ ground tours separately with `python -m scheduling.tplus rat99-m5 --output-dir routes_tplus`, replacing `rat99-m5` with the desired source route name. |
| §5.6 | Run `python run_suite.py endurance endurance --plans local_runs/main/full` after generating the main solutions. The runner replays the matching locally generated 15% reference schedule when one is available. |
| §5.7 | Run `python run_suite.py main release --plans local_runs/main/full`. The runner reports cases without a matching reference schedule. Results depend on the supplied schedule, including its parcel batches and flight sequence. See the [perturbation protocol](datasets/release_uncertainty/protocol.json). |
| §5.8 | Run `python compare_management.py --runs local_runs/main/full --out local_runs/management` to compare paired objective bounds and the marginal benefit of each additional resupply drone. |

The four published manifests contain 216 main, 24 truck-only, 432 F/T makespan, and 864 endurance input records. The makespan command covers the published F/T inputs. The T+ utility produces ground tours, and the paper's T+ and sufficient-capacity makespans require additional scheduling runs. Local time-limited runs and regenerated reference schedules can yield different results from the paper's archived runs.

## Inputs and documentation

The fixed truck routes and delivery-drone sorties come from the public route archive accompanying Ren, Froger, Jabali, and Liang (2024). The fixed truck-only comparison uses six CVRPLIB route sets. The [source-route catalog](datasets/source_routes.json) records coordinates, routes, sorties, and source identifiers. [Ren et al. (2024)](https://doi.org/10.1016/j.ejor.2024.05.031) · [CVRPLIB](https://galgos.inf.puc-rio.br/cvrplib/en/instances)

See the [input guide](datasets/README.md) for formats and units, the [algorithm guide](ALGORITHMS.md) for solver options, and the [historical benchmarks](datasets/legacy/README.md) for earlier benchmarks and extraction commands.
