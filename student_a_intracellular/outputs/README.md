# outputs: one folder per task

| Folder | Task | Main files |
|---|---|---|
| `a1/` | A1: CRP iModulons from PRECISE-1K (Lamoureux 2023 Fig. 4) | `crp_activity_by_condition.csv` (Crp-1 / Crp-2 activity per sample), `condition_response_curve.png` |
| `a2/` | A2: causal-mediation test of CRP → gene edges | `mediation_results.csv`, `mediation_summary.json`, `a2_mediation_report.md`, `fig_a2_edges.png` |
| `a3/` | A3: CRP/lac circuit model, fit and hysteresis check | `circuit_model.yaml` (+ SBML `circuit_model.xml`), `validated_edges.json` (every circuit edge and its evidence), `parameter_fit_report.md`, `hysteresis_check.png`, `fig_a3_*.png` |
| `a4/` | A4: stochastic simulator, Module 1 (Vilar 2003) | `a4_report.md`, `summary_stats.json`, trajectories (`*.parquet`), `fig_a4_*.png` |
| `a5/` | A5: noise calibration on PETRI-seq | `calibrated_circuit_model.yaml` (the model A6 and Student B use), `calibration_report.md`, `fig_*.png` |
| `a6/` | A6: perturbation interface, Module 2, and DiMA check (Lamoureux 2023 Fig. 6) | `perturbation_dima_comparison.csv`, `a6_report.md`, `fig_a6_*.png`, `heldout/` |

Each task reads the folder before it: A3 reads `a1/` and `a2/`, A4 reads `a3/circuit_model.yaml`, A5 writes `a5/calibrated_circuit_model.yaml` from it, and A6 and Student B load that calibrated model.
