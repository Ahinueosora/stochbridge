# student_a_intracellular — data-grounded intracellular engine (StochBridge Track A)

E. coli CRP / lac circuit, from bulk transcriptomes (PRECISE-1K) and single-cell RNA-seq (PETRI-seq) to a calibrated stochastic
single-cell simulator (Module 1) with a perturbation interface (Module 2). Student B's code reads `outputs/` and imports
`src/ssa_simulator.py` and `src/perturbation.py`; the formats are documented in `../shared/schemas/`.

## Layout
```
data/precise1k/      PRECISE-1K + K-12 modulome files (not in git; see "Data" below)
data/scrnaseq/GSE141018/   PETRI-seq summary extracted from GEO (Blattman et al. 2020)
src/                 code (one entry point per task, listed below)
config/              crp_lac_circuit.yaml (A3 human-readable spec), precise1k_condition_map.yaml (A6)
notebooks/           01 … 05, one per task, executed
outputs/             what Student B reads; per-task folders a4/, a5/, a6/
tests/               pytest, fixed seeds (38 tests)
```

## Tasks, entry points and outputs

| Task | What it does | Run | Main outputs |
|---|---|---|---|
| A1 iModulons | Reproduce Crp-1 / Crp-2 iModulons (Lamoureux 2023 Fig. 4) from PRECISE-1K; CRP activity across conditions | `notebooks/01_crp_reproduction.ipynb` | `outputs/condition_response_curve.png/.csv` |
| A2 causal edges | Causal-mediation test of CRP → target edges (gbmMINER Fig. 5B–D logic), checked against RegulonDB | `python src/mediation_test.py` (`notebooks/a2_ecoli_miner_causal_pipeline.ipynb` is the earlier exploratory MINER analysis) | `outputs/a2/mediation_results.csv`, `mediation_summary.json`, `outputs/validated_edges.json` |
| A3 circuit_model | Lac core (Ozbudak 2004) + Kuhlman 2007 CRP terms + 21 CRP outputs; fit and hysteresis | `python src/run_a3.py`, `python src/build_a3_v3.py` | `outputs/circuit_model.yaml/.xml`, `parameter_fit_report.md`, `hysteresis_check.png`; `notebooks/02_hysteresis_check.ipynb` |
| A4 stochastic_engine (Module 1) | Gillespie SSA (GillesPy2), Vilar 2003 all-or-none induction | `python src/a4_vilar.py` | `outputs/a4/` (trajectories, summary_stats.json, report); `notebooks/03_vilar_reproduction.ipynb` |
| A5 noise_calibration | Fit extrinsic noise to PETRI-seq, validate exp/stationary split | `python src/load_scrnaseq.py data/scrnaseq/GSE141018/GSM4489500_*.txt.gz` (rebuilds the summary), then `python src/fit_noise_params.py && python src/a5_validate.py && python src/a5_figures.py && python src/a5_sensitivity.py` | `outputs/calibrated_circuit_model.yaml`, `outputs/calibration_report.md`; `notebooks/04_noise_calibration.ipynb` |
| A6 perturbation_interface (Module 2) | Perturbation classes; Lamoureux 2023 Fig. 6 DiMA vs matched real conditions | `python src/a6_reproduce_aat.py && python src/a6_dima_compare.py && python src/a6_heldout.py && python src/a6_figures.py` | `outputs/perturbation_dima_comparison.csv`, `outputs/a6/`; `notebooks/05_dima_workflow.ipynb` |

## Module 1 — stochastic_engine (`src/ssa_simulator.py`, `src/io_schema.py`)
```python
from ssa_simulator import StochasticCircuit
from perturbation import InducerLevelChange
sc = StochasticCircuit("outputs/calibrated_circuit_model.yaml")      # compiles once
traj = sc.run(InducerLevelChange(30, n_cells=1000, seed=1))          # -> Trajectory (long table: cell_id, time, *_count)
traj.to_parquet("traj.parquet"); traj.write_summary("summary_stats.json")
```
Modes: `fast` (internal inducer at quasi-steady state; default) and `exact` (all reactions). Schema: `../shared/schemas/trajectory_schema.md`.

## Module 2 — perturbation_interface (`src/perturbation.py`)
```python
from perturbation import load_perturbations, GeneKnockout, InducerLevelChange, CarbonSourceShift
modified = GeneKnockout("crp").apply(model_dict)                    # new model dict, input untouched
for pert in load_perturbations("config/precise1k_condition_map.yaml"):
    traj = StochasticCircuit(pert.apply(model_dict)).run(seed=1)
```
Contract: `../shared/schemas/perturbation_schema.md`.

## Data
Copy from the PRECISE-1K repository (`SBRG/precise1k`, `data/`) into `data/precise1k/` under these names:
`p1k_precise1k.json.gz`, `p1k_log_tpm_qc.csv`, `p1k_metadata_qc.csv`, `p1k_counts.csv` (from `precise1k/`),
`k12_counts.csv`, `k12_metadata_qc.csv` (from `k12_modulome/`).
PETRI-seq: GEO GSE141018 (GSM4489500); the raw count matrix and the summary used by A5 are in `data/scrnaseq/GSE141018/`.

## Tests
`pytest -q tests` — 38 tests: mediation rule and the 21 edges (A2), circuit checks (A3), fixed-seed SSA regression and fast-vs-exact (A4), noise layer (A5),
perturbation API and Fig. 6 reproduction (A6).
