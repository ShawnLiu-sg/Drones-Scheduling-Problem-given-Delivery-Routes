# Algorithm entry points

Use Python 3.11 with Gurobi 12.0.1 and a working Gurobi license. Install the packages in requirements.txt. The interval, full threshold, makespan and endurance procedures can invoke IBM ILOG CP Optimizer (12.10 in the study). Set the environment variable CP_OPTIMIZER to the executable path or make cpoptimizer available on PATH. SciPy and NumPy are used for release uncertainty.

Run a single instance as follows, replacing the two paths with the problem JSON and a local output file.

    python run_solver.py full INPUT.json --output local/full.json
    python run_solver.py liu INPUT.json --output local/liu.json
    python run_solver.py release INPUT.json --plan local/full.json --output local/risk.json

The solver accepts the input JSON, reads the truck stops, delivery sorties and parcel releases, and constructs the no-resupply baseline and a deterministic starting plan. An explicit input H is a hard horizon. With no hard horizon, the event-network size and maximum lag define a safe horizon. A feasible total-delay incumbent can further tighten the delay optimization horizon. Time budgets default to 120 seconds for small and medium cases and 300 for large cases. Use --seconds to override. Gurobi and CP Optimizer use one worker. Outputs include the input identity and input hash. The JSON plan includes each truck's arrival and departure times, delivery launches and recoveries, and parcel assignments and replenishment dispatch times.

| Method | Model or algorithm |
| --- | --- |
| full | Optional-interval primal search, A3 interval-capacity bounds with four-tick delay segments and dynamic observation intervals, then cumulative threshold search |
| base | Threshold search with both optional components disabled |
| no_interval | Only the interval-capacity strengthening component |
| no_capacity | Only optional-interval search |
| liu | Adapted depot network flow. Empty sortie lists select pure-truck delivery |
| compact, interval | Cumulative threshold formulation or optional replenishment interval formulation in isolation |
| graph_a | Scheme A chain-selection master with difference-graph Benders cuts |
| price_d | Scheme D branch-and-price |
| price_original, price_filter, price_signature, price_closure, price_all, price_core, price_adaptive | Original, domain filtering, timing signature, sparse closure, combined and projected pricing variants |
| price_milp | Capacity-strengthened master with MILP pricing |
| fixed_policy | Time-indexed reference formulation with a prescribed parcel-to-stop assignment |
| reference, pilot_milp | Time-indexed reference formulation and compact pilot MILP |
| cbs, cbs_dp | Fleet conflict search with the local route oracle or dynamic-programming oracle |
| makespan | Absolute threshold algorithm for latest truck or shared-drone return |
| endurance | Fixed-plan replay, retiming, new dispatch construction, physical feasibility and absolute-return-threshold feasibility |
| release | Monte Carlo and analytical probabilities under four laws and three perturbation scales, given a locally generated plan |

In the capacity stage, the solver first computes the layered load bound. If that bound leaves the global gap open without improving the entry lower bound, it uses the remaining capacity-stage time to search for common-stop batches. A completed local proof can strengthen an interval coefficient or exclude a delay segment. A local search that reaches its time or node limit retains the layered bound. The result records the gate decision under `stats.energy_details.batch_gate` when the capacity stage runs.

For a historical prescribed-handoff comparison, extract the policy and run the fixed-policy model:

    python datasets/legacy/extract_case.py synthetic_controlled controlled-k3-n5-s001-j2-a06-m2 policy_case.json --fixed-policy
    python run_solver.py fixed_policy policy_case.json --output local/policy.json

The model fixes every parcel's handoff station and optimizes the dispatch and synchronization times. It constructs its starting plan using those prescribed handoffs. Only a physically verified plan satisfying the policy can tighten the search horizon. For a paired free-handoff solve, extract the same series and label without the extraction flag and run the desired free-handoff method.

For experiment 4 run makespan separately on each arm. Use --plan with a compatible locally generated plan to supply an additional upper bound. The release package contains inputs only. This local warm-start procedure may produce different limited-time witnesses from an archived run.

For experiment 5 run endurance on an endurance configuration. With --plan, supply a plan generated on the matching 15% instance. The replay freezes that plan before permitting rescheduling. The reference-return target defaults to the sum of bare truck returns, identical across the formal endurance arms. Use --reference-total to set a different absolute target. The returned physical and operational fields carry separate conclusions. For experiments 4 and 5, pass --evidence-dir NEW_LOCAL_DIR to retain the generated model and solver logs, or omit it to store them temporarily.

For experiment 6 the input issue scenario, parcel order, releases, weights and fixed plan must match. A plan can be produced with full or liu and consumes no external result file. Its physical replay precedes stochastic simulation. Each arm reports the simulated probability, Wilson interval and independent product-of-CDF probability. Recomputing a reference plan can change the conditional probability.

Experiment 7 uses the bounds from local optimizations. The functions pair_status and triple_bounds in scheduling/comparisons.py compute the signed differences and marginal-gain certificates from endpoint lower bounds, upper bounds and statuses. The alpha, late-parcel fraction and fleet comparisons must match the original route and other parameters. Additional components include deterministic truck-route insertion, positive-cycle certificate construction, a warm-start improvement heuristic and a small two-team exhaustive oracle.
