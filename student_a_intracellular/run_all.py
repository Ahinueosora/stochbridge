"""
Re-run every task A1 -> A6 from the raw data, in order, and re-execute the notebooks.

    python get_data.py          # once: downloads PRECISE-1K (pinned version)
    python run_all.py           # everything, about 25 min on a 2-core machine (+40 min if outputs/a4/exact_*.parquet are deleted)
    python run_all.py a2 a3     # only some tasks (each task reads the outputs of the one before)
    python verify.py            # then check the results against the numbers in the report

Every step writes into outputs/aN/ and overwrites the files that are in the repository, so after a full run
`git status` / `git diff` shows exactly what changed. Simulations use fixed seeds.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable
DATA = "data/precise1k/precise1k/data"


def nb(name):
    return [PY, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace",
            "--ExecutePreprocessor.timeout=3600", f"notebooks/{name}"]


STEPS = {
    "a1": [nb("01_crp_reproduction.ipynb")],
    "a2": [[PY, "src/mediation_test.py"]],
    "a3": [[PY, "src/run_a3.py", DATA, "outputs/a3"],
           [PY, "../../src/build_network.py", "../../data/curated/lac_network_gene_check.xlsx", {"cwd": "outputs/a3"}],
           [PY, "src/build_a3_v3.py", "outputs/a3/network_v3_edges.csv", "outputs/a3"],
           nb("02_hysteresis_check.ipynb")],
    "a4": [[PY, "src/a4_vilar.py"], [PY, "src/a4_validate_modes.py", "400"], nb("03_vilar_reproduction.ipynb")],
    "a5": [[PY, "src/fit_noise_params.py"], [PY, "src/a5_validate.py"], [PY, "src/a5_figures.py"],
           [PY, "src/a5_sensitivity.py"], nb("04_noise_calibration.ipynb")],
    "a6": [[PY, "src/a6_reproduce_aat.py"], [PY, "src/a6_dima_compare.py"], [PY, "src/a6_heldout.py"],
           [PY, "src/a6_figures.py"], nb("05_dima_workflow.ipynb")],
}


def main(tasks):
    if not (ROOT / DATA / "precise1k" / "log_tpm_qc.csv").exists():
        sys.exit("Data missing: run `python get_data.py` first.")
    # the notebook kernel must be this same Python (with the pinned packages), even if the venv is not activated
    env = dict(os.environ, MPLBACKEND="Agg", PATH=os.path.dirname(PY) + os.pathsep + os.environ.get("PATH", ""))
    t_all = time.time()
    for task in tasks:
        for step in STEPS[task]:
            opts = step[-1] if isinstance(step[-1], dict) else {}
            cmd = [s for s in step if not isinstance(s, dict)]
            cwd = ROOT / opts.get("cwd", ".")
            shown = " ".join(c if c != PY else "python" for c in cmd)
            print(f"[{task}] {shown}", flush=True)
            t0 = time.time()
            r = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if r.returncode != 0:
                print(r.stdout[-4000:])
                sys.exit(f"FAILED at [{task}] {shown}")
            print(f"      done in {time.time() - t0:.0f} s", flush=True)
    print(f"All requested tasks finished in {(time.time() - t_all) / 60:.1f} min. Now run: python verify.py")


if __name__ == "__main__":
    wanted = [a.lower() for a in sys.argv[1:]] or list(STEPS)
    bad = [t for t in wanted if t not in STEPS]
    if bad:
        sys.exit(f"Unknown task(s) {bad}; choose from {list(STEPS)}")
    main(wanted)
