# Module 2 · perturbation_interface (Student A, A6) — contract for Student B's screen_runner (B3)

Import path: `student_a_intracellular/src/perturbation.py` (put `src/` on `sys.path`).

```python
import yaml
from perturbation import load_perturbations, Perturbation, GeneKnockout, InducerLevelChange, CarbonSourceShift, compose
from ssa_simulator import StochasticCircuit

model = yaml.safe_load(open("outputs/calibrated_circuit_model.yaml"))
perts = load_perturbations("config/precise1k_condition_map.yaml")            # list[Perturbation], file order
# with_controls=True -> [(perturbation, in_silico_control, matched_real_condition_or_None), ...]

for pert in perts:
    modified = pert.apply(model)          # dict, deep copy; `model` is untouched
    sc = StochasticCircuit(modified)      # compiles once per model (~seconds)
    for rep in range(1, n_reps + 1):
        traj = sc.run(seed=rep)           # runs the embedded perturbation; any run field can be overridden here
        traj.to_parquet(f"{pert.name}_rep{rep}.parquet")
```

Faster equivalent (one compile for the whole screen, identical trajectories): `sc = StochasticCircuit(model); sc.run(pert, seed=rep)`.

## Perturbation fields (all serialisable: `to_dict()` / `from_dict()` / `to_json()`)

| field | meaning | default |
|---|---|---|
| name, kind, description | identifiers (`kind` = class name, or "A+B" for composed) | "baseline" |
| TMG_uM, IPTG_uM | inducer outside the cell (TMG via LacY; IPTG passive) | 0 |
| glucose_uM, carbon_source | glucose seen by the circuit; carbon source label | 0, "" |
| deleted_genes, crp_knockout | knockouts (crp/cyaA set crp_knockout) | [], False |
| param_overrides | any circuit parameter, e.g. {"rho": 5} | {} |
| initial_state | "uninduced" or "induced" history | "uninduced" |
| t_end_min, dt_min, n_cells, seed | run settings (seed must be >= 1) | 1200, 5, 100, 1 |

## What `apply()` writes into the model
`parameters` (overrides), `genetics` {deleted_genes, crp_knockout}, `perturbation` (the full dict above), `name` suffixed
`__<perturbation name>`. Allowed knockouts: crp, cyaA, lacI, lacY, and the 21 CRP-output genes; anything else raises
`NotImplementedError` (it would have no effect inside the circuit). PTS sugars other than glucose raise `NotImplementedError`.

## Validation hook
Every entry in `config/precise1k_condition_map.yaml` with `matched:` names a real PRECISE-1K / K-12 condition
(`project:condition`). `src/a6_dima_compare.py` turns simulated populations into pseudo-bulk samples and runs the
Lamoureux 2023 Fig. 6 DiMA against it; results in `outputs/perturbation_dima_comparison.csv`.
