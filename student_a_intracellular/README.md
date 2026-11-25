# student_a_intracellular: data-grounded intracellular engine (StochBridge Track A)

*E. coli* CRP / *lac* circuit, built from bulk transcriptomes (PRECISE-1K) and single-cell RNA-seq (PETRI-seq), and turned
into a calibrated stochastic single-cell simulator (Module 1) with a perturbation interface (Module 2).
Student B's code reads `outputs/` and imports `src/ssa_simulator.py` and `src/perturbation.py`; the formats are in `../shared/schemas/`.

Anuradha Barnwal (2023BB10659), BBD451 B.Tech Major Project, IIT Delhi. Supervisor: Prof. Ishaan Gupta.

---

## Verify the results (about 10 minutes)

Needs **Python 3.11 or newer** (tested with 3.13 on Linux) and **git**. The simulator (GillesPy2) compiles a small C++ program, so a
C++ compiler is needed: on macOS run `xcode-select --install` once; on Ubuntu/Debian `sudo apt install g++`.

```bash
git clone https://github.com/Ahinueosora/stochbridge
cd stochbridge/student_a_intracellular
python3 -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt                            # exact versions that produced the results
python get_data.py                                         # downloads PRECISE-1K (~240 MB, pinned commit, checksums checked)
python verify.py                                           # recomputes, tests and checks every headline number
```

`verify.py` prints one line per check and ends with a count. Expected ending:

```
  PASS  A6  held-out refit (CRP curves refitted without the 32 compared samples): 9 of 11 correct vs 10 of 12  [held out 9/11, in-sample 10/12]
25 of 25 checks passed.
```
(`python verify.py --a3` adds three A3 checks: `28 of 28 checks passed.`) Any `FAIL` line says what differed.

What `verify.py` does:

1. **Data**: the 10 PRECISE-1K files match the SHA-256 checksums of the pinned version; the PETRI-seq file is present.
2. **Rerun from raw data**: A1 (notebook 01) and A2 (causal-mediation test) are recomputed in a temporary folder and compared
   with `outputs/`: they must be identical. `python verify.py --a3` also recomputes A3 (+3 min).
3. **Tests**: the 38 pytest tests, including fixed-seed simulator regressions and the Lamoureux Fig. 6 reproduction.
4. **Results**: each task's headline numbers in `outputs/` are compared with the numbers quoted in the report
   (for example: 21 causal edges of 543 genes; switch-off thresholds within 23%; switch-on at 37.3 µM;
   Lamoureux Fig. 6 found 11/11; 16 of 22 simulated DiMA calls correct).

### Regenerate everything from scratch

```bash
python run_all.py          # all tasks A1 -> A6 and the notebooks, about 25 min on a 2-core machine
python run_all.py a2 a3    # or only some tasks; each task reads the outputs of the one before
python verify.py
git status                 # shows which output files changed (with the pinned versions, none should change in content)
```

With the pinned package versions every step is deterministic (all simulations use fixed seeds), so `run_all.py`
reproduces the files in `outputs/`; this was checked by regenerating every output from a fresh clone. Two kinds of file
are expected to differ: figure PNGs (matplotlib embeds metadata) and wall-clock timings stored in the JSON files.
The exact-mode SSA runs used for the A4 fast-vs-exact check take about 40 minutes, so they are kept in
`outputs/a4/exact_*.parquet` and reused; delete those three files to recompute them (the result is identical). Newer package versions also run, but can change numbers slightly:
for example pandas 3.1 changed how missing values turn into text, which is why the code now handles missing values explicitly.

---

## Layout
```
get_data.py          downloads the PRECISE-1K files (pinned commit 4829b83 of SBRG/precise1k, SHA-256 checked)
run_all.py           reruns every task in order
verify.py            checks data, reruns A1-A2 (A3 with --a3), runs the tests, compares headline numbers with the report
data/precise1k/      PRECISE-1K + K-12 modulome files (downloaded by get_data.py, not in git)
data/scrnaseq/       PETRI-seq count matrix from GEO GSE141018 (GSM4489500), in git (3.9 MB)
data/curated/        hand-curated inputs, in git: RegulonDB evidence per CRP edge, the combined network table
src/                 code (entry points per task below)
config/              crp_lac_circuit.yaml (A3 human-readable spec), precise1k_condition_map.yaml (A6)
notebooks/           01 ... 05, one per task, executed
outputs/             one folder per task: a1/ ... a6/, each with a README explaining every figure and file (index: outputs/README.md)
tests/               pytest, fixed seeds (38 tests)
```

## Tasks, entry points and outputs

`run_all.py` runs exactly these commands (from this folder).

| Task | What it does | Commands | Time* | Main outputs |
|---|---|---|---|---|
| A1 iModulons | Reproduce Crp-1 / Crp-2 iModulons (Lamoureux 2023 Fig. 4) from PRECISE-1K; CRP activity per sample | `jupyter nbconvert --execute --inplace notebooks/01_crp_reproduction.ipynb` | 10 s | `outputs/a1/` |
| A2 causal edges | Causal-mediation test of CRP → target edges (gbmMINER Fig. 5B–D logic) | `python src/mediation_test.py` | 1 min | `outputs/a2/` (`mediation_results.csv`, `mediation_summary.json`, `a2_mediation_report.md`) |
| A3 circuit model | *lac* core (Ozbudak 2004) + Kuhlman 2007 CRP terms + 21 CRP outputs; fit and hysteresis | `python src/run_a3.py data/precise1k/precise1k/data outputs/a3`<br>`(cd outputs/a3 && python ../../src/build_network.py ../../data/curated/lac_network_gene_check.xlsx)`<br>`python src/build_a3_v3.py outputs/a3/network_v3_edges.csv outputs/a3` | 3 min | `outputs/a3/` (`circuit_model.yaml/.xml`, `validated_edges.json`, `parameter_fit_report.md`, `hysteresis_check.png`); notebook 02 |
| A4 Module 1 | Gillespie SSA (GillesPy2); Vilar 2003 all-or-none induction; fast vs exact mode | `python src/a4_vilar.py`<br>`python src/a4_validate_modes.py 400` | 14 min (+40 min if the cached exact runs `outputs/a4/exact_*.parquet` are deleted) | `outputs/a4/` (`a4_report.md`, `vilar_results.json`, `summary_stats.json`, trajectories); notebook 03 |
| A5 noise calibration | Fit cell-to-cell noise to PETRI-seq; validate the growth-phase split | `python src/fit_noise_params.py`<br>`python src/a5_validate.py`<br>`python src/a5_figures.py`<br>`python src/a5_sensitivity.py` | 5 min | `outputs/a5/` (`calibrated_circuit_model.yaml`, `calibration_report.md`); notebook 04 |
| A6 Module 2 | Perturbation classes; Lamoureux 2023 Fig. 6 DiMA vs matched real conditions | `python src/a6_reproduce_aat.py`<br>`python src/a6_dima_compare.py`<br>`python src/a6_heldout.py`<br>`python src/a6_figures.py` | 3 min | `outputs/a6/` (`perturbation_dima_comparison.csv`, `a6_report.md`, `heldout/`); notebook 05 |

\*On a 2-core cloud machine.

`notebooks/a2_ecoli_miner_causal_pipeline.ipynb` is the earlier exploratory MINER analysis (it needs the separate MINER
package and is not part of the pipeline above). `python src/load_scrnaseq.py data/scrnaseq/GSE141018/GSM4489500_*.txt.gz`
rebuilds the PETRI-seq summary (`petri_2_01_summary_regenerated.json`) from the raw GEO matrix; a test checks it matches the one A5 uses.

## Module 1: stochastic_engine (`src/ssa_simulator.py`, `src/io_schema.py`)
```python
# run from src/ (or add src/ to PYTHONPATH)
import yaml
from ssa_simulator import StochasticCircuit
from perturbation import GeneKnockout, InducerLevelChange

model = yaml.safe_load(open("../outputs/a5/calibrated_circuit_model.yaml"))
sc = StochasticCircuit(model)                                            # compiles once
traj = sc.run(InducerLevelChange(30, n_cells=1000, seed=1))             # -> Trajectory (cell_id, time, *_count)
traj.to_parquet("traj.parquet"); traj.write_summary("summary_stats.json")
ko = StochasticCircuit(GeneKnockout("crp").apply(model)).run(seed=1)    # crp deletion
```
Modes: `fast` (internal inducer at quasi-steady state; default) and `exact` (all reactions). Schema: `../shared/schemas/trajectory_schema.md`.

## Module 2: perturbation_interface (`src/perturbation.py`)
```python
from perturbation import load_perturbations, GeneKnockout, InducerLevelChange, CarbonSourceShift
modified = GeneKnockout("crp").apply(model)                       # new model dict, input untouched
for pert in load_perturbations("../config/precise1k_condition_map.yaml"):
    traj = StochasticCircuit(pert.apply(model)).run(seed=1)
```
Contract: `../shared/schemas/perturbation_schema.md`.

## Data sources
- **PRECISE-1K and the public K-12 modulome**: Lamoureux et al. 2023, *Nucleic Acids Res* 51:10176, github.com/SBRG/precise1k
  (commit 4829b83, 29 Aug 2023). Downloaded by `get_data.py`.
- **PETRI-seq**: Blattman et al. 2020, *Nat Microbiol* 5:1192; GEO GSE141018, sample GSM4489500. In `data/scrnaseq/`.
- **RegulonDB** release 14.5.1 (accessed 24 Nov 2026) and the combined network table: in `data/curated/` (see its README).

## Tests
`python -m pytest -q tests`: 38 tests covering the mediation rule and the 21 edges (A2), circuit checks (A3), fixed-seed SSA
regression and fast vs exact (A4), the noise layer (A5), and the perturbation API and Fig. 6 reproduction (A6).
