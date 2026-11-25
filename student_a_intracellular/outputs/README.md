# outputs: one folder per task

Every folder has a **README.md** that explains each figure (what it shows, how to read it, the outcome) and lists the files,
a **figures/** subfolder with the key figures, and the result files written by the task's code.

| Folder | Task | Paper reproduced | Outcome | Main result files |
|---|---|---|---|---|
| [`a1/`](a1) | A1: CRP iModulons from PRECISE-1K | Lamoureux et al. 2023, Fig. 4 | **Strong.** Reproduced exactly (R² = 1.00, all 1,035 samples) | `crp_activity_by_condition.csv`, `condition_response_curve.png` |
| [`a2/`](a2) | A2: causal-mediation test of CRP → gene edges | Turkarslan et al. 2024 (gbmMINER), Fig. 5B–D | **Good, small sample.** 21 causal edges of 543 genes; 17 agree with RegulonDB | `mediation_results.csv`, `mediation_summary.json`, `a2_mediation_report.md` |
| [`a3/`](a3) | A3: CRP/lac circuit model, fit and hysteresis | Ozbudak et al. 2004 | **Strong.** Held-out switch-off thresholds within ~23%; cusp ρ = 9.01 | `circuit_model.yaml` (+ SBML `.xml`), `validated_edges.json`, `parameter_fit_report.md` |
| [`a4/`](a4) | A4: stochastic simulator (Module 1) | Vilar et al. 2003 | **Strong.** All-or-none, bimodal; switch-less strain unimodal | `vilar_results.json`, `summary_stats.json`, `fast_vs_exact.json`, `a4_report.md` |
| [`a5/`](a5) | A5: noise calibration on PETRI-seq | Blattman et al. 2020 | **Adequate.** Distributions match; false switching 19% → 0.9%; phase split weaker | `calibrated_circuit_model.yaml` (used by A6 and Student B), `calibration_report.md` |
| [`a6/`](a6) | A6: perturbation interface (Module 2), DiMA check | Lamoureux et al. 2023, Fig. 6 | **Adequate.** Fig. 6 reproduced 11/11; 16 of 22 simulated calls correct; effects ~3.6× too small | `perturbation_dima_comparison.csv`, `perturbation_summary.csv`, `a6_report.md`, `heldout/` |

Each task reads the folder before it: A3 reads `a1/` and `a2/`, A4 reads `a3/circuit_model.yaml`, A5 writes
`a5/calibrated_circuit_model.yaml` from it, and A6 and Student B load that calibrated model.
`python ../run_all.py` regenerates all of them; `python ../verify.py` checks them.
