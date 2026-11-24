# Scheduling Replenishment Drones on Fixed Delivery Routes

This repository contains algorithms and benchmark instances for scheduling a shared fleet of depot-based replenishment drones. Each truck follows a fixed route. Its onboard delivery drone serves customers on predetermined sorties and rejoins the truck. When a parcel becomes available after the truck has left the depot, a replenishment drone carries it to a stop on the assigned truck's route. The scheduler chooses which parcels travel together, where to hand them over, and when to dispatch each replenishment flight so that truck and drone operations stay synchronized.

## Experiments

| Experiment | Instances | What to run |
| --- | --- | --- |
| 1. Total return delay | 216 scenarios with fixed truck routes and delivery sorties | Compare `full` and `liu` on the same inputs. The objective is to minimize the sum of truck return delays relative to the baseline routes without replenishment. |
| 2. Algorithm ablation | 120 main scenarios with fewer replenishment drones than trucks | Compare `base`, `no_capacity`, `no_interval`, and `full` on the same subset. |
| 3. Truck-only delivery routes | 24 CVRPLIB scenarios | Compare `full` and `liu` when trucks serve all customers and replenishment drones supply the trucks. |
| 4. Joint makespan | 216 matched pairs, 432 configurations | Run `makespan` on the truck-and-drone and truck-only delivery plans for the same customers. Completion occurs when all trucks and all dispatched replenishment drones have returned to the depot. |
| 5. Onboard delivery drone endurance | 216 scenarios at endurance margins of 0%, 5%, 15%, and 30%, for 864 configurations | Use `endurance` to replay a local 15% reference schedule, then allow rescheduling under the selected endurance margin. |
| 6. Parcel release uncertainty | Main scenarios with a locally verified schedule, four error distributions, and three perturbation levels | Use `release` to estimate the probability that the fixed schedule remains feasible and compare the estimate with its analytical value. |
| 7. Release timing and replenishment capacity | Matched main scenarios that vary one factor at a time, plus triples with one, two, and three replenishment drones | Run `compare_management.py` on locally computed objective bounds to compare delay changes and the marginal benefit of each additional drone. |

The four formal input manifests contain 1,536 configurations. Scenarios shared across experiments appear in each applicable manifest. The repository also includes earlier graph and pricing algorithms (schemes A and D), threshold variants, conflict search, small-instance enumeration, and historical benchmarks. See the [algorithm guide](ALGORITHMS.md), [formal input guide](datasets/README.md), and [historical benchmark guide](datasets/legacy/README.md).

## Run locally

Run the commands below from the repository root. Use Python 3.11 and a licensed installation of Gurobi. Some methods also require IBM ILOG CP Optimizer. Add `cpoptimizer` to your `PATH` or set `CP_OPTIMIZER` to the path of its executable.

To validate the inputs, compare two solvers on one small instance, and evaluate the resulting schedule under release uncertainty, run these commands.

```sh
python -m pip install -r requirements.txt
python datasets/validate_inputs.py
python run_suite.py main full --tier small --first 1
python run_suite.py main liu --tier small --first 1
python run_suite.py main release --tier small --first 1 --plans local_runs/main/full
```

`run_suite.py` processes each manifest sequentially and writes outputs to `local_runs/<dataset>/<method>/`. Omit `--first` to run every instance in the selected size tier, and omit `--tier` to cover all sizes. The default time limit per instance is 120 seconds for small and medium instances and 300 seconds for large instances. Use `--ablation` to select the experiment 2 subset. For single-instance methods and options, run `python run_solver.py --help`.

The following commands run the main suite, the truck-only comparison method, the paired makespan study, and the analyses that use locally generated schedules or bounds.

```sh
python run_suite.py main full
python run_suite.py cvrp liu
python run_suite.py makespan makespan
python run_suite.py endurance endurance --plans local_runs/main/full
python run_suite.py main release --plans local_runs/main/full
python compare_management.py --runs local_runs/main/full --out local_runs/management
```

Experiment 6 evaluates the particular schedule supplied with `--plans`. Two schedules with the same objective value can have different replenishment departure times or parcel batches, so their probabilities of remaining feasible can differ. The batch runner lists scenarios that lack a local reference schedule.

## Data sources

The fixed truck routes and onboard delivery sorties for the main experiments come from the public archive accompanying Ren, Froger, Jabali, and Liang (2024). The truck-only comparison uses six route sets from CVRPLIB. The route input file records node coordinates, stop sequences, delivery sorties, and source identifiers. [Route inputs](datasets/source_routes.json) · [Ren et al. (2024)](https://doi.org/10.1016/j.ejor.2024.05.031) · [CVRPLIB](https://galgos.inf.puc-rio.br/cvrplib/en/instances)
