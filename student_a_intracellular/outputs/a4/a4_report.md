# A4 report: stochastic simulator (Module 1) and the Vilar 2003 reproduction

**Result: the SSA of our own circuit (A3 v3) gives all-or-none induction.** Single cells switch from OFF to ON suddenly and at random times; the population splits into two separate groups (bimodal, dip test p < 0.001), while a graded strain with the same average stays one group (p > 0.9). Cells that start ON keep their memory; cells that start OFF drift ON slowly inside the memory zone, which Vilar's model did not show. That drift is the main thing A5 must calibrate.

## 1. What was built

| File | Role |
|---|---|
| `src/ssa_simulator.py` | `StochasticCircuit(yaml)` loads `outputs/a3/circuit_model.yaml` (A3 v3), compiles once (GillesPy2 C++ SSA, `variable=True`), `run(perturbation) -> Trajectory` |
| `src/io_schema.py` | the Trajectory contract: columns, types, metadata, `validate`, parquet I/O, summary stats (mean, variance, Fano, ON fraction, bimodality) |
| `src/perturbation.py` | `Perturbation` + the four brief classes: `InducerLevelChange`, `CarbonSourceShift`, `GeneKnockout('crp')`, `OperatorDecoy(rho)` (A6 extends) |
| `shared/schemas/trajectory_schema.md` | the contract for Student B (version 1.0.0, change policy) |
| `tests/test_ssa_reproducibility.py` | 12 tests: fixed seeds bit-identical, golden regression, schema, and behaviour (basal = ODE, memory, bimodal vs unimodal, Δcrp OFF) |
| `notebooks/03_vilar_reproduction.ipynb` | the reproduction, executed (reads cached runs; one line reruns everything) |
| `src/a4_vilar.py`, `src/a4_validate_modes.py` | simulations behind the notebook |
| `outputs/a4/summary_stats.json` | canonical Module 1 summary (the 37.3 µM run, 2,000 cells, 30 h) |
| `outputs/a4/trajectory_*.parquet` | example Trajectories (2,000 cells × 361 time points incl. 21 CRP-target mRNAs) |

**Two modes, same output columns.**
- `exact`: all 8 reactions, internal TMG simulated molecule by molecule (~4.5 s per cell-20 h when switching, ~30 s when ON).
- `fast` (default): internal TMG held at its quasi-steady level inside the repression reaction. It relaxes in ~2 min while LacY changes over hours, and it numbers in the thousands. 1,000 cells × 20 h in 1–30 s.

**Fast vs exact (independent seeds; 120 cells each; 30 for cells starting ON):**

| Case | % ON at 20 h, exact | % ON at 20 h, fast | Switch-time distributions |
|---|---|---|---|
| 40 µM, cells start OFF | 47.5 | 52.5 (diff −4 ± 6 pp) | KS p = 0.97 |
| 15 µM, cells start OFF | 15.0 | 10.8 (diff +4 ± 5 pp) | KS p = 0.36 |
| 15 µM, cells start ON | 100 | 100; ON LacY 1,399 vs 1,403, Fano 45 vs 42 | – |

No detectable difference, so the population runs use fast mode; the molecular-level figure uses exact mode.

## 2. Vilar et al. 2003 Fig. 3, with our circuit

Doses scaled to our memory zone (ODE 4.1–30 µM; Vilar used 7 µM sub-optimal, 5 µM maintenance, 500 µM high in his model).

| Vilar panel | Ours | Result |
|---|---|---|
| 3a sub-optimal dose: cells switch abruptly at random times; smooth average; population curve below it | 37.3 µM, 2,000 cells, 30 h | One cell goes 10% → 90% ON in **300 min** (median); the 2,000-cell average takes **1,256 min**. Switch times spread from 380 to 1,165 min (IQR). 67.5% switched by 30 h. Population curve (ON cells grow 15% slower, assumed) ends at 557 vs average 937 LacY |
| 3b maintenance dose: no switching from either history | 11 µM, 1,000 cells per history, 30 h | Cells starting ON: **0 / 1,000 flip OFF** (as Vilar). Cells starting OFF: **17.4% flip ON** within 30 h (Vilar: none) |
| 3c high dose: everyone switches; average ≈ population | 500 µM, 2,000 cells | 99.4% ON by 30 h, median switch 90 min; average 1,492 vs population 1,439 |
| molecular level | exact SSA, one cell at 40 µM | repressor falls off → mRNA bursts → LacY → more TMG inside → promoter stays free (`fig_a4_molecular.png`) |

**Bimodal, not just noisy** (`fig_a4_bimodality.png`): at 37.3 µM the LacY histogram has two peaks from 8 h on (dip test p < 0.001 at 8, 12, 20 and 30 h; Sarle coefficient 0.75–0.89). At 4 h only 10% of cells are ON and the test is not yet significant (p = 0.91). An `OperatorDecoy(rho = 5)` strain at a dose that gives the same 30 h average (0.18 µM, mean 899 vs 937) has one peak at every time (dip p 0.77–0.99; Sarle 0.38–0.56). The dip test is run on de-quantised counts; without that, the spikes at 0, 1, 2 molecules read as extra peaks.

## 3. SSA hysteresis (1,000 cells per dose and history, 20 h)

| | SSA, 50% of cells | SSA, 95% | ODE, 20 h | Analytic | Ozbudak |
|---|---|---|---|---|---|
| Switch ON (cells start OFF) | **37 µM** | not reached by 200 µM (91%) | 36 µM | 30 µM | 30 µM |
| Switch OFF (cells start ON) | **3.4 µM** | 2.2 µM | 2.9 µM | 4.1 µM | 3 µM |

- Coming from ON, the switch is sharp: 24% → 78% → 98% ON between 2.9 and 5.9 µM.
- Coming from OFF, it is broad: 7% at 12 µM, 33% at 24 µM, 46% at 34 µM, 91% at 200 µM. Cells with very little LacY must wait for a burst before uptake can start, and the model has no passive TMG entry (Vilar includes one).

## 4. Findings for A5 and the panel

1. **OFF cells leak ON inside the memory zone** (7% at 12 µM in 20 h; 17% at 11 µM in 30 h). This is burst-driven: one mRNA makes ~20 LacY (k_tl / d_m, a set value), enough to start the feedback near the upper threshold. The LacI fall-off rate `k_unrep` barely matters (6.5–7.6% across 0.02–0.5 per min with ρ fixed). A5 should calibrate **burst size** and `x0` against PETRI-seq / single-cell data, not `k_unrep`.
2. **Switching from OFF is slow and broad**; adding passive TMG entry (Vilar's b₂ term) is the obvious model change to test.
3. **Δcrp stays OFF** at 500 µM TMG in the SSA (mean LacY < 2), as in Kuhlman and Traxler.
4. LacY Fano factor in fully ON cells is ~45 (both modes); the A5 noise target.

## 5. Caveats

- Cells do not divide: dilution is a continuous loss (as in Vilar's cellular level). Division noise would add to the ON-state Fano factor.
- The population curve uses an assumed growth cost of 15% for ON cells (Vilar used measured rates).
- Units of LacY come from A3's chosen rates (~1,500 when fully ON); A5 recalibrates.
- The golden regression values depend on GillesPy2 1.8.3 and its C++ random-number generator; regenerate deliberately with `python tests/test_ssa_reproducibility.py --update-golden`.

## 6. Reproduce

```bash
pip install gillespy2 pyarrow pyyaml scipy pandas matplotlib diptest python-libsbml
python src/a4_validate_modes.py 120      # fast vs exact (~25 min, exact mode is slow)
python src/a4_vilar.py                   # all Vilar / hysteresis runs (~15 min on 2 cores)
python -m pytest -q tests                # 23 tests (A3 + A4)
jupyter nbconvert --execute --to notebook --inplace notebooks/03_vilar_reproduction.ipynb
```

## Update 24 Nov 2026: larger fast-vs-exact check (400 cells; 100 for cells starting ON)

| Case | Switched by 20 h, exact | Switched by 20 h, fast | Comparison |
|---|---|---|---|
| 40 µM, cells start OFF | 51.8% | 55.0% | difference −3.3 ± 3.5 points; switch times KS p = 0.87 |
| 15 µM, cells start OFF | 14.3% | 18.8% | difference −4.5 ± 2.6 points; switch times KS p = 0.56 |
| 15 µM, cells start ON | 100% stay ON | 100% stay ON | ON LacY 1,415 vs 1,400; Fano 35 vs 40 |

- **What agrees.** Switch times and the levels and noise of ON cells match between the two modes.
- **Small bias.** In both cases that start OFF, the fast mode switches slightly more cells ON. Pooled difference: 4.1 ± 2.1 points (p ≈ 0.05). The 120-cell check could not resolve it.
- **What it means.** Switching fractions from the fast mode may be a few points high.
- **Files.** Results in `outputs/a4/fast_vs_exact.json`; the 120-cell version is kept as `fast_vs_exact_n120.json`.
- **Bug found and fixed.** This run exposed a bug from A6: a perturbation subclass could not be rebuilt from its own `to_dict()` output. The classes now round-trip, and a test checks it.
