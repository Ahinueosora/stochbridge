"""
A5 robustness checks.

1. Which genes set eta2: all 16 fitted operons vs only the 15 operons the model simulates (circuit outputs + lacZYA),
   and leave-one-gene-out (jackknife).
2. The assumed capture per molecule q (0.05; range 0.025-0.10): eta2 does not depend on q (Eq. for B); the burst bound
   and the validation do, so the validation is rerun at q = 0.025 and 0.10.

Writes outputs/a5/sensitivity.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import fit_noise_params as fnp  # noqa: E402

ROOT = HERE.parent
SIMULATED = {"lacZYA", "hflC", "tdcB", "uhpT", "cas3", "rbs", "yahN", "fecE", "mhpF", "flhDC", "araJ", "nadC", "galK",
             "gyrA", "malK", "preA", "idnO"}


def fit_subset(S, keep, cv2k):
    S2 = json.loads(json.dumps(S))
    S2["operons"] = {k: v for k, v in S["operons"].items() if k in keep or k.startswith("_")}
    f = fnp.fit_fano_line(S2, "exp", cv2k)
    return {"n_genes": f["n_genes"], "eta2": f["eta2"], "eta2_se": f["eta2_se"], "genes": sorted(f["genes"])}


def main():
    S = fnp.load_summary()
    cap = json.loads((ROOT / "outputs" / "a5" / "noise_fit.json").read_text())["capture"]["exp"]["cv2"]
    base = fnp.fit_fano_line(S, "exp", cap)
    used = set(base["genes"])
    R = {"all_fitted_operons": {"n_genes": base["n_genes"], "eta2": base["eta2"], "eta2_se": base["eta2_se"], "genes": sorted(used)},
         "simulated_operons_only": fit_subset(S, used & SIMULATED, cap)}
    jack = {g: fit_subset(S, used - {g}, cap)["eta2"] for g in sorted(used)}
    R["leave_one_gene_out"] = {"min": min(jack.values()), "max": max(jack.values()), "by_gene_left_out": jack}
    R["burst_bound_by_q"] = base["burst_b_by_q0"]
    # validation at other capture rates
    from a5_validate import run_model
    val = {}
    for q in (0.025, 0.05, 0.10):
        m = run_model(ROOT / "outputs" / "calibrated_circuit_model.yaml", S, f"q={q}", q0=q)
        a = m["agreement"] if "agreement" in m else None
        from a5_validate import compare
        ag, _ = compare(S, m)
        val[str(q)] = {"pca17_auc": m["pca17"]["auc"], "exp_hist_chi2_per_dof": ag["exp"]["hist_chi2"] / ag["exp"]["hist_dof"],
                       "exp_fano_rmse": ag["exp"].get("fano_rmse_expressed"),
                       "stat_median_abs_log2_mean_ratio": ag["stat"]["median_abs_log2_mean_ratio"]}
    R["validation_by_q"] = val
    (ROOT / "outputs" / "a5" / "sensitivity.json").write_text(json.dumps(R, indent=2, default=float))
    return R


if __name__ == "__main__":
    R = main()
    print(json.dumps({k: (v if k != "leave_one_gene_out" else {kk: vv for kk, vv in v.items() if kk != "by_gene_left_out"})
                      for k, v in R.items() if k != "burst_bound_by_q"}, indent=1, default=float)[:3000])
