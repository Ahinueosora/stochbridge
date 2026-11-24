"""
A3 add-on: the 21 A2-validated CRP edges as extra outputs of the circuit.

For each target, expression is a Hill function of CRP-cAMP activity c, fitted on PRECISE-1K:
    activated:  TPM = b + V * c^h / (K^h + c^h)
    repressed:  TPM = b + V * K^h / (K^h + c^h)
Direction is taken from A2 (knockout lowers target -> activated).

CRP activity c is read from the gat operon, the CRP-only promoter Ozbudak used as his CRP
reporter: c = gat TPM / median gat TPM of MG1655 non-glucose samples (so c = 1 is the
no-glucose state of the A3 circuit, and c = c_crp(G) links directly to the lac model).

Validation:
  1. In PRECISE-1K: predicted knockout effect (crp project) vs A2's measured total effect.
  2. Out of sample: K-12 modulome 'raman' project, cyaA knockout + cAMP 0/0.1/0.5/1 mM + WT,
     a direct intervention on CRP-cAMP that the fit never saw.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

# A2 causal edges: name -> total effect (log2, knockout/deletion vs WT)
A2_EDGES = {
    "hflC": 0.45, "tdcB": -0.97, "uhpT": -2.37, "cas3": -0.67, "rbsC": -1.94, "yahN": -1.32,
    "fecE": 1.03, "mhpF": 0.81, "flhC": -3.00, "flhD": -3.12, "araJ": 0.53, "rbsA": -2.53,
    "nadC": 0.32, "galK": -0.69, "rbsB": -2.34, "gyrA": 0.34, "rbsR": -1.20, "lamB": -4.48,
    "malK": -4.71, "preA": -0.84, "idnO": -0.87,
}
A2_KO_EFFECT = {  # knockout-only column from the A2 table
    "hflC": 0.59, "tdcB": -1.22, "uhpT": -2.89, "cas3": -0.85, "rbsC": -2.86, "yahN": -1.69,
    "fecE": 1.57, "mhpF": 0.87, "flhC": -3.91, "flhD": -4.04, "araJ": 0.78, "rbsA": -3.92,
    "nadC": 0.34, "galK": -0.81, "rbsB": -3.63, "gyrA": 0.44, "rbsR": -1.75, "lamB": -5.88,
    "malK": -6.17, "preA": -1.19, "idnO": -1.32,
}
GAT = ["b2096", "b2095", "b2094"]  # gatY, gatZ, gatA
TPM_TO_MOLECULES = 0.0024           # ~2,400 mRNAs per cell (assumption, for stochastic export)


def hill(c, b, V, K, h, sign):
    c = np.maximum(np.asarray(c, float), 1e-6)
    act = c**h / (K**h + c**h)
    return b + V * (act if sign > 0 else 1 - act)


def crp_meter_p1k(root: Path):
    d = root / "precise1k"
    md = pd.read_csv(d / "metadata_qc.csv", index_col=0)
    lt = pd.read_csv(d / "log_tpm_qc.csv", index_col=0)
    keep = md.Strain != "DGF-298"                      # genome-reduced strain, lacks many genes
    md, lt = md[keep], lt.loc[:, keep]
    gat = (2 ** lt.loc[GAT] - 1).mean()
    cs = md["Carbon Source (g/L)"].astype(str)
    ref_mask = (md.Strain == "MG1655") & ~cs.str.startswith("glucose") & (cs != "nan") \
        & ~md["Strain Description"].astype(str).str.contains("crp|cya", case=False)
    ref = gat[ref_mask].median()
    return md, lt, gat / ref, ref


def fit_targets(root: Path, outdir: Path):
    gi = pd.read_csv(root / "annotation" / "gene_info.csv", index_col=0)
    name2b = {n: b for b, n in gi.gene_name.items() if n in A2_EDGES}
    md, lt, c, ref = crp_meter_p1k(root)
    rows, curves = [], {}
    for name, eff in A2_EDGES.items():
        b_id = name2b[name]
        y = lt.loc[b_id].values                       # log2(TPM+1)
        sign = -1 if eff > 0 else +1                  # KO raises target -> CRP represses it
        tpm = 2 ** y - 1

        def resid(th):
            b, V, lK, h = th
            return np.log2(hill(c.values, b, V, 10 ** lK, h, sign) + 1) - y

        x0 = [np.percentile(tpm, 10), np.percentile(tpm, 90), 0.0, 1.0]
        sol = optimize.least_squares(resid, x0, bounds=([0, 0, -3, 0.5], [1e5, 1e6, 2, 4]))
        b, V, lK, h = sol.x
        r2 = 1 - np.sum(sol.fun ** 2) / np.sum((y - y.mean()) ** 2)
        rho_s = stats.spearmanr(c.values, y).correlation
        # prediction of the A2 knockout effect: crp project, KO vs WT
        crp = md.project == "crp"
        ko = crp & md["Strain Description"].str.contains("del_crp")
        wt = crp & ~md["Strain Description"].str.contains("crp_|del_crp")
        pred = lambda m: np.log2(hill(c[m].values, b, V, 10 ** lK, h, sign) + 1).mean()  # noqa: E731
        rows.append({"target": name, "b_number": b_id, "direction": "activated" if sign > 0 else "repressed",
                     "basal_TPM": b, "amplitude_TPM": V, "K_crp": 10 ** lK, "hill_n": h,
                     "R2_all_samples": r2, "spearman_vs_crp": rho_s,
                     "a2_ko_effect_log2": A2_KO_EFFECT[name], "pred_ko_effect_log2": pred(ko) - pred(wt)})
        curves[name] = (b, V, 10 ** lK, h, sign)
    tab = pd.DataFrame(rows)
    tab.to_csv(outdir / "crp_outputs_fit.csv", index=False)
    return tab, curves, ref


def validate_raman(root: Path, curves, ref, tpm_k12: pd.DataFrame | None = None):
    """Out-of-sample: cyaA KO + cAMP titration (K-12 modulome)."""
    md = pd.read_csv(root / "k12_modulome" / "metadata_qc.csv", index_col=0)
    s = md[md.project == "raman"]
    if tpm_k12 is None:
        cnt = pd.read_csv(root / "k12_modulome" / "counts.csv", index_col=0)[s.index]
        gi = pd.read_csv(root / "annotation" / "gene_info.csv", index_col=0)
        L = (gi.end - gi.start).abs().reindex(cnt.index) / 1000
        ok = L.notna() & (L > 0)
        rpk = cnt[ok].div(L[ok], axis=0)
        tpm_k12 = rpk / rpk.sum() * 1e6
    else:
        tpm_k12 = tpm_k12[s.index]
    gi = pd.read_csv(root / "annotation" / "gene_info.csv", index_col=0)
    name2b = {n: b for b, n in gi.gene_name.items() if n in curves}
    c = tpm_k12.loc[GAT].mean() / ref
    grp = s.condition
    rows = []
    for name, (b, V, K, h, sign) in curves.items():
        obs = np.log2(tpm_k12.loc[name2b[name]] + 1)
        prd = np.log2(hill(c.values, b, V, K, h, sign) + 1)
        df = pd.DataFrame({"obs": obs.values, "pred": prd, "cond": grp.values})
        m = df.groupby("cond").mean()
        for cond in m.index:
            if cond == "wt":
                continue
            rows.append({"target": name, "condition": cond,
                         "obs_log2FC_vs_wt": m.loc[cond, "obs"] - m.loc["wt", "obs"],
                         "pred_log2FC_vs_wt": m.loc[cond, "pred"] - m.loc["wt", "pred"]})
    v = pd.DataFrame(rows)
    return v, c.groupby(grp).median()


def run(root: Path, outdir: Path, tpm_k12=None):
    outdir.mkdir(parents=True, exist_ok=True)
    tab, curves, ref = fit_targets(root, outdir)
    v, c_by = validate_raman(root, curves, ref, tpm_k12)
    v.to_csv(outdir / "crp_outputs_validation_raman.csv", index=False)
    r_ko = stats.pearsonr(tab.a2_ko_effect_log2, tab.pred_ko_effect_log2)
    r_ra = stats.pearsonr(v.obs_log2FC_vs_wt, v.pred_log2FC_vs_wt)
    sign_ra = float(np.mean(np.sign(v.obs_log2FC_vs_wt) == np.sign(v.pred_log2FC_vs_wt)))
    summary = {
        "gat_reference_TPM": float(ref),
        "median_R2": float(tab.R2_all_samples.median()),
        "ko_effect_pearson_r": float(r_ko[0]), "ko_effect_p": float(r_ko[1]),
        "ko_effect_slope_pred_over_obs": float(np.polyfit(tab.a2_ko_effect_log2, tab.pred_ko_effect_log2, 1)[0]),
        "raman_pearson_r": float(r_ra[0]), "raman_p": float(r_ra[1]), "raman_sign_agreement": sign_ra,
        "raman_crp_activity_by_condition": c_by.round(3).to_dict(),
    }
    (outdir / "crp_outputs_summary.json").write_text(json.dumps(summary, indent=2))
    return tab, v, summary, curves


if __name__ == "__main__":
    tab, v, summary, _ = run(Path(sys.argv[1]), Path(sys.argv[2]))
    pd.set_option("display.width", 200)
    print(tab.round(3).to_string())
    print(json.dumps(summary, indent=1))
