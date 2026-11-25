"""
A6 held-out test for the PRECISE-1K comparisons.

The A3 CRP-output curves were fitted on all PRECISE-1K samples, which include the strains A6 compares against
(projects `crp` and `cra_crp`, and the reference wt_glc samples). Here the 21 curves are refitted with those samples
removed, written into a copy of the calibrated model, and the four PRECISE-1K perturbations are rerun. Any agreement
that survives is then out of sample, like the K-12 (cyaA, IPTG) comparisons.

Outputs: outputs/a6/heldout/{crp_outputs_fit_heldout.csv, calibrated_circuit_model_heldout.yaml,
         perturbation_dima_comparison.csv, perturbation_summary.csv, ...}
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import optimize

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from crp_outputs import A2_EDGES, GAT, hill  # noqa: E402
from dima import P1K_DIR, P1K_REFS, KnowledgeBase  # noqa: E402

ROOT = HERE.parent
OUTD = ROOT / "outputs" / "a6" / "heldout"
HELD_PROJECTS = {"crp", "cra_crp"}
P1K_PERTURBATIONS = ["glucose_to_glycerol", "glucose_to_acetate", "dcrp_glucose", "dcrp_glycerol"]


def refit(exclude):
    kb = KnowledgeBase.load(with_log_tpm=False)
    md = pd.read_csv(P1K_DIR / "p1k_metadata_qc.csv", index_col=0)
    lt = pd.read_csv(P1K_DIR / "p1k_log_tpm_qc.csv", index_col=0)
    keep = (md.Strain != "DGF-298") & ~md.index.isin(exclude)
    md, lt = md[keep], lt.loc[:, keep]
    gat = (2 ** lt.loc[GAT] - 1).mean()
    cs = md["Carbon Source (g/L)"].fillna("nan").astype(str)
    ref_mask = (md.Strain == "MG1655") & ~cs.str.startswith("glucose") & (cs != "nan") \
        & ~md["Strain Description"].fillna("nan").astype(str).str.contains("crp|cya", case=False)
    c = gat / gat[ref_mask].median()
    name2b = {n: b for b, n in kb.gene_table.gene_name.items() if n in A2_EDGES}
    rows = []
    for name, eff in A2_EDGES.items():
        y = lt.loc[name2b[name]].values
        sign = -1 if eff > 0 else +1
        tpm = 2 ** y - 1

        def resid(th):
            b, V, lK, h = th
            return np.log2(hill(c.values, b, V, 10 ** lK, h, sign) + 1) - y

        sol = optimize.least_squares(resid, [np.percentile(tpm, 10), np.percentile(tpm, 90), 0.0, 1.0],
                                     bounds=([0, 0, -3, 0.5], [1e5, 1e6, 2, 4]))
        b, V, lK, h = sol.x
        rows.append({"target": name, "basal_TPM": b, "amplitude_TPM": V, "K_crp": 10 ** lK, "hill_n": h,
                     "R2": 1 - np.sum(sol.fun ** 2) / np.sum((y - y.mean()) ** 2)})
    return pd.DataFrame(rows).set_index("target"), int(keep.sum()), int(ref_mask.sum())


def main():
    OUTD.mkdir(parents=True, exist_ok=True)
    md = pd.read_csv(P1K_DIR / "p1k_metadata_qc.csv", index_col=0)
    exclude = set(md.index[md.project.isin(HELD_PROJECTS)]) | set(P1K_REFS) \
        | set(md.index[(md.project == "control") & (md.condition == "wt_glc")])
    tab, n_fit, n_ref = refit(exclude)
    tab.to_csv(OUTD / "crp_outputs_fit_heldout.csv")
    model = yaml.safe_load((ROOT / "outputs" / "a5" / "calibrated_circuit_model.yaml").read_text())
    old = copy.deepcopy(model["crp_outputs"]["targets"])
    for name, r in tab.iterrows():
        t = model["crp_outputs"]["targets"][name]
        t.update({"basal_TPM": float(r.basal_TPM), "amplitude_TPM": float(r.amplitude_TPM), "K_crp": float(r.K_crp),
                  "hill_n": float(r.hill_n)})
    model["name"] = model["name"] + "_heldout_crp_fit"
    model["crp_outputs"]["note_heldout"] = (f"CRP-output curves refitted without the {len(exclude)} samples A6 compares "
                                            f"against (projects {sorted(HELD_PROJECTS)} and control:wt_glc)")
    path = OUTD / "calibrated_circuit_model_heldout.yaml"
    path.write_text(yaml.safe_dump(model, sort_keys=False, width=130))
    from a6_dima_compare import main as compare
    df, S = compare(model=path, a6_dir=OUTD, main_csv=OUTD / "perturbation_dima_comparison.csv", only=P1K_PERTURBATIONS)
    shift = {n: float(np.log2((tab.loc[n].amplitude_TPM + 1) / (old[n]["amplitude_TPM"] + 1))) for n in tab.index}
    info = {"excluded_samples": len(exclude), "samples_used_for_fit": n_fit, "reference_samples_for_c": n_ref,
            "median_R2_heldout_fit": float(tab.R2.median()),
            "log2_amplitude_change_vs_full_fit": shift}
    (OUTD / "heldout_info.json").write_text(json.dumps(info, indent=2))
    return S, info


if __name__ == "__main__":
    S, info = main()
    print(json.dumps({k: v for k, v in info.items() if k != "log2_amplitude_change_vs_full_fit"}, indent=1))
    cols = ["perturbation", "gene_r", "gene_sign_agreement_moved", "panel_real_dimas", "panel_sim_dimas", "panel_recovered",
            "panel_fig6_real_dimas", "panel_fig6_recovered", "panel_multigene_real_dimas", "panel_multigene_recovered", "panel_r",
            "r_sim_vs_ceiling"]
    print(S[cols].round(2).to_string())
