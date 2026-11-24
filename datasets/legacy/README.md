# Historical scheduling benchmarks

The index maps 2,256 experimental labels from 20 benchmark series to 2,107 distinct physical inputs. Repeated uses of the same physical input share one JSON line. These series include the initial A/D benchmarks, controlled and pressure families, Solomon-derived scenarios, public-route pilots, the nested fleet and package experiments, 0.5 kg and 0.1 kg payload variants, and finite-horizon verification cases. Route arrays describe the fixed truck and delivery-drone operations supplied to the scheduling problem. The current seven experiments are listed separately in the main datasets.

Each row of `instances.jsonl` contains an identifier, scenario descriptors, a load unit, and the solver's physical instance. The first materialization of a shared instance supplies those descriptors. Select the original experiment by the pair of columns `series` and `original_label` in `index.csv`. The `physical_id` joins the experiment label to its input.

From the repository root, extract one case into a new local file:

```shell
python datasets/legacy/extract_case.py initial_A_D k2-p3-m1-s0 case.json
```

- `M` counts the available shared replenishment drones. The payload capacity `Q` and parcel weights `w` use the same units. Multiply those values by `weight_unit_kg` where a physical conversion is supplied. A null conversion denotes abstract benchmark load units.
- `travel`, `service`, `tau` give truck arc durations, station service durations and one-way depot flight durations. For each delivery sortie, `l` and `b` are launch and recovery stop indices, and `F`, `E` bound the elapsed off-truck time. Parcel `dest` selects a truck stop or a sortie, depending on its `kind`. Release times `r` are absolute. Some checks intentionally contain incompatible parameters.
- An explicit `instance.H` occurs in historical finite-horizon validation cases; the horizon constrains the problem as posed. Main benchmark instances use the route and resource data. Their old computational search windows came from algorithm-specific incumbents or construction bounds. Recompute a valid window before any time-indexed solve. This can change the computational search domain and timing relative to the original runs.
- The historical fixed-policy comparisons prescribed the demand station or a fixed alternating sequence of stops before scheduling. The applicable rule is marked per experiment row. Use `--fixed-policy` when extracting that restricted variant. The default extracts the unrestricted scheduling instance.
- Under `solomon_derived`, fixed-route travel data were derived from C101, R101 and RC101 from the [SINTEF Solomon collection](https://www.sintef.no/projectweb/top/vrptw/100-customers/). The archived public-route pilots used the published route families identified in the scenario descriptors, including the Ren, Schmidt and Xu collections. Synthetic controls and hand-built cases use their original benchmark identifiers.
