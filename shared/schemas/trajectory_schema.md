# Trajectory schema (Module 1 → Module 4)

**Owner:** Student A (stochastic_engine, A4). **Consumer:** Student B, Module 4 (per agent, per time step).
**Version:** 1.0.0 (24 Nov 2026). **Code:** `student_a_intracellular/src/io_schema.py` (`Trajectory`, `validate`, `SCHEMA_VERSION`).

> **Change policy.** The core columns, their order, names, types and units are frozen for all 1.x versions.
> Adding an optional column is a minor bump (1.0 → 1.1). Renaming, removing, retyping or changing units is a major
> bump (→ 2.0) and **must be agreed with Student B before it is merged**.

## What one file is

One `.parquet` file = one run of one `Perturbation`: `n_cells` independent stochastic cells (one cell = one agent),
all sampled on the same time grid. Long format: one row per (cell, time).

## Columns

| # | Column | Type | Unit | Meaning |
|---|---|---|---|---|
| 1 | `cell_id` | int32 | – | agent index, 0 … n_cells − 1 |
| 2 | `time` | float64 | min | time since the perturbation started; grid 0, dt, 2·dt, …, t_end (same for every cell) |
| 3 | `P_free_count` | int64 | molecules | lac promoter free of LacI (0 or 1) |
| 4 | `P_rep_count` | int64 | molecules | lac promoter bound by LacI (0 or 1); `P_free + P_rep = 1` |
| 5 | `M_count` | int64 | molecules | lacZYA mRNA |
| 6 | `Y_count` | int64 | molecules | LacY permease: **the lac switch read-out** |
| 7 | `X_count` | int64 | molecules | internal TMG (simulated, or derived in fast mode; see `x_column`) |
| 8+ | `mRNA_<gene>_count` | int64 | molecules | optional: the 21 A2-validated CRP targets (hflC, tdcB, …), present only if the run used `include_outputs=True` |

Rows are sorted by `cell_id`, then `time`. Counts are never negative. `io_schema.validate(df)` enforces all of this.

## File metadata (parquet key-value, JSON strings)

| Key | Content |
|---|---|
| `schema_version` | `"1.0.0"` |
| `units` | `{"time": "min", "counts": "molecules per cell"}` |
| `model` | name, version and SHA-256 of the `circuit_model.yaml` used |
| `perturbation` | the full Perturbation (TMG_uM, glucose_uM, crp_knockout, param_overrides, initial_state, t_end_min, dt_min, n_cells, seed, kind, name) |
| `solver`, `gillespy2_version` | GillesPy2 SSACSolver (Gillespie direct method, C++) |
| `mode` | `"exact"` (all reactions) or `"fast"` (internal TMG at quasi-steady state) |
| `x_column` | `"simulated"` or `"qssa_derived"` (X = k_up_eff · Y / d_x) |
| `seed` | integer; same seed + same file versions ⇒ identical file |
| `rates`, `initial_counts` | the numbers actually used |
| `Y_on`, `Y_off`, `on_threshold_Y` | deterministic ON and OFF LacY for this condition; a cell is **ON** if `Y_count > on_threshold_Y` (geometric mean of the two) |

## Companion summary (`summary_stats*.json`)

Per species, across cells, at the last time point: `mean`, `variance`, `fano`; the same over the second half of the
run; `fraction_ON`; bimodality of log10(LacY + 1) (Sarle coefficient, Hartigan dip test).

## How Module 4 reads it

```python
from io_schema import Trajectory
tr = Trajectory.from_parquet("trajectory_vilar3a_sub.parquet")   # validates on load
tr.at(600)                     # all agents at t = 600 min
tr.wide("Y")                   # agents x time matrix of LacY
on = tr.at(600).Y_count > tr.meta["on_threshold_Y"]
```

Plain pandas/pyarrow also works: `pd.read_parquet(path)`; metadata via `pyarrow.parquet.read_schema(path).metadata`.

## A5 addition: noise metadata
Trajectories from `calibrated_circuit_model.yaml` carry `meta["noise"]` (output_extrinsic_cv2, output_burst_b, crp_extrinsic_sigma,
lac_protein_burst) and `meta["crp_outputs_method"]` ("exact immigration-death sampling, per-cell gamma rate" when η² > 0, else "SSA").
`mRNA_<gene>_count` columns are molecules, not UMIs; to compare with scRNA-seq apply `observation_model` (Poisson thinning with per-cell capture).
