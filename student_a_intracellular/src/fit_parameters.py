"""
A3 - parameter fitting for the CRP/lac circuit.

Two data sources, each constraining what it can measure:
  1. PRECISE-1K + K-12 modulome RNA-seq  -> repression factor rho, CRP scaling of lac,
     glucose effect on CRP activity.  (Bulk data: no hysteresis, no bimodality.)
  2. Ozbudak et al. 2004 Fig. 2c (digitized) -> uptake scale and glucose half-effects,
     fitted to SWITCH-ON thresholds only; SWITCH-OFF thresholds are held out as the test.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).parent))
import reaction_network as rn  # noqa: E402

GENES = {"lacZ": "b0344", "lacY": "b0343", "lacA": "b0342", "lacI": "b0345",
         "gatY": "b2096", "gatZ": "b2095", "gatA": "b2094"}

# Ozbudak Fig 2c, digitized by eye from the PDF (log axes; read error ~ +/-15%).
# G = 0 values are from the text (30 and 3 uM).
FIG2C = pd.DataFrame({
    "glucose_uM": [0, 10, 25, 50, 100, 200, 1000],
    "switch_on_uM": [30, 165, 400, 580, 650, 690, 790],
    "switch_off_uM": [3, 17, 43, 65, 94, 107, 140],
})


# =============================================================== 1. RNA-seq
def load_precise1k(root: Path):
    d = root / "precise1k"
    md = pd.read_csv(d / "metadata_qc.csv", index_col=0)
    lt = pd.read_csv(d / "log_tpm_qc.csv", index_col=0)
    A = pd.read_csv(d / "A.csv", index_col=0)
    df = pd.DataFrame({k: 2 ** lt.loc[v] - 1 for k, v in GENES.items()})   # log2(TPM+1) -> TPM
    df["gat"] = df[["gatY", "gatZ", "gatA"]].mean(axis=1)
    df["Crp1"], df["Crp2"] = A.loc["Crp-1"], A.loc["Crp-2"]
    df = df.join(md[["project", "condition", "Strain", "Strain Description", "Carbon Source (g/L)", "Evolved Sample"]])
    sd = df["Strain Description"].astype(str)
    wt = ~sd.str.contains("crp|cya|lac|gat|ptsG|ptsH|ptsI|crr", case=False)
    wt &= df["Evolved Sample"].astype(str).str.lower() != "yes"
    # Lac-intact background only: BW25113 carries lacZ4787 deletion (+ lacI^q) and DGF-298 has
    # lost the whole lac region; both bias uninduced lac downward. Other strains (GMOS, W3110)
    # are lac-intact but heterogeneous, so the reference is MG1655 only.
    wt &= df["Strain"] == "MG1655"
    return df[wt].copy()


def load_k12_srna(root: Path):
    d = root / "k12_modulome"
    md = pd.read_csv(d / "metadata_qc.csv", index_col=0)
    s = md[md.project == "srna_oe"]
    c = pd.read_csv(d / "counts.csv", index_col=0)
    c = c[[x for x in s.index if x in c.columns]]
    g = pd.read_csv(root / "annotation" / "gene_info.csv", index_col=0)
    L = (g.end - g.start).abs().reindex(c.index) / 1000
    ok = L.notna() & (L > 0)
    rpk = c[ok].div(L[ok], axis=0)
    tpm = rpk / rpk.sum() * 1e6
    out = pd.DataFrame({k: tpm.loc[v] for k, v in GENES.items()})
    out["gat"] = out[["gatY", "gatZ", "gatA"]].mean(axis=1)
    out = out.join(s[["condition", "Strain Description", "Supplement"]])
    out["iptg"] = out.Supplement.astype(str).str.contains("IPTG")
    return out


def rho_from_rnaseq(p1k, k12, n_boot=5000, seed=1):
    """rho = (lac/gat)_induced / (lac/gat)_uninduced.  gat is the CRP-only promoter Ozbudak
    used as a CRP meter (HcRed); dividing by it removes CRP differences between samples."""
    rng = np.random.default_rng(seed)
    induced = k12[k12.iptg]
    unind_batch = k12[k12.condition == "wt"]
    unind_p1k = p1k[p1k.gat > 50]
    rows = []
    for gene in ["lacY", "lacZ"]:
        ri = (induced[gene] / induced.gat).values
        for ref, u in [("same batch (K-12 srna_oe WT, n=3)", unind_batch),
                       (f"PRECISE-1K MG1655 WT (n={len(unind_p1k)})", unind_p1k)]:
            ru = (u[gene] / u.gat).values
            est = np.median(ri) / np.median(ru)
            bs = [np.median(rng.choice(ri, len(ri))) / np.median(rng.choice(ru, len(ru)))
                  for _ in range(n_boot)]
            lo, hi = np.percentile(bs, [2.5, 97.5])
            rows.append({"gene": gene, "uninduced_reference": ref, "rho": est, "ci_lo": lo, "ci_hi": hi})
    per_group = (induced.assign(r=induced.lacY / induced.gat).groupby("condition").r.median())
    return pd.DataFrame(rows), per_group


def crp_scaling(p1k):
    """Glucose effect on CRP activity and on basal lac (PRECISE-1K, WT-like strains)."""
    cs = p1k["Carbon Source (g/L)"].astype(str)
    glc, other = p1k[cs.str.startswith("glucose")], p1k[~cs.str.startswith("glucose") & (cs != "nan")]
    res = {}
    for k in ["gat", "lacY", "lacZ"]:
        res[k] = {"median_glucose_TPM": float(glc[k].median()), "median_other_TPM": float(other[k].median()),
                  "fold_other_over_glucose": float(other[k].median() / glc[k].median()),
                  "mwu_p": float(stats.mannwhitneyu(glc[k], other[k]).pvalue),
                  "n_glucose": int(len(glc)), "n_other": int(len(other))}
    m = (p1k.lacY > 0.5) & (p1k.gat > 1)
    r = stats.linregress(np.log2(p1k.loc[m, "gat"]), np.log2(p1k.loc[m, "lacY"]))
    res["basal_lacY_vs_gat_loglog"] = {"slope": r.slope, "se": r.stderr, "r": r.rvalue, "n": int(m.sum())}
    return res


# =============================================================== 2. Fig 2c
def calibrate_k_up(p, T_on0=30.0):
    """Set uptake scale so the wild-type switch-on threshold at 0 glucose is 30 uM."""
    p = dict(p, crp_K=p.get("crp_K") or 1.0, excl_K=p.get("excl_K") or 1.0, k_up=1.0)
    _, a_on = rn.switching_drives(p["rho"], p["n"])
    p["k_up"] = a_on / rn.drive(T_on0, 0.0, p)
    return p["k_up"]


def fit_glucose(p, fit_rho=False):
    """Fit glucose half-effects (and optionally rho) to Fig 2c.
    Default: switch-on thresholds only (switch-off held out)."""
    data = FIG2C[FIG2C.glucose_uM > 0]
    targets = ["switch_on_uM"] + (["switch_off_uM"] if fit_rho else [])

    def unpack(theta):
        q = dict(p, crp_K=10 ** theta[0], excl_K=10 ** theta[1])
        if fit_rho:
            q["rho"] = 10 ** theta[2]
        q["k_up"] = calibrate_k_up(q)
        return q

    def resid(theta):
        q = unpack(theta)
        if rn.switching_drives(q["rho"], q["n"]) is None:
            return np.full(len(data) * len(targets), 5.0)
        out = []
        for _, row in data.iterrows():
            th = rn.thresholds_tmg(row.glucose_uM, q)
            if th is None:                  # CRP-weakened rho fell below the cusp at this glucose
                out += [5.0] * len(targets)
                continue
            off, on = th
            pred = {"switch_on_uM": on, "switch_off_uM": off}
            out += [np.log10(min(pred[t], 1e6)) - np.log10(row[t]) for t in targets]
        return np.array(out)

    x0 = [1.5, 1.0] + ([np.log10(p["rho"])] if fit_rho else [])
    sol = optimize.least_squares(resid, x0, bounds=([-1, -1] + ([0.96] if fit_rho else []),
                                                     [4, 4] + ([4] if fit_rho else [])))
    q = unpack(sol.x)
    return q, sol


def threshold_table(p):
    rows = []
    for _, r in FIG2C.iterrows():
        off, on = rn.thresholds_tmg(r.glucose_uM, p) or (np.nan, np.nan)
        rows.append({"glucose_uM": r.glucose_uM, "on_obs": r.switch_on_uM, "on_model": on,
                     "off_obs": r.switch_off_uM, "off_model": off, "rho_eff": rn.rho_eff(r.glucose_uM, p)})
    t = pd.DataFrame(rows)
    t["on_log10_err"] = np.log10(t.on_model / t.on_obs)
    t["off_log10_err"] = np.log10(t.off_model / t.off_obs)
    return t


def run_all(data_root: Path, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    p1k = load_precise1k(data_root)
    k12 = load_k12_srna(data_root)
    rho_tab, per_group = rho_from_rnaseq(p1k, k12)
    crp = crp_scaling(p1k)

    # Primary model: rho = 170 (Ozbudak, GFP protein; inside the RNA-seq envelope)
    p = dict(rn.DEFAULT_PARAMS)
    p, sol = fit_glucose(p)
    tab = threshold_table(p)
    # Sensitivity: rho at the RNA-seq estimates
    sens = {}
    for label, extra in [("rho_66_batch_matched_lacY", {"rho": 66.0}),
                         ("rho_435_PRECISE1K_MG1655_lacZ", {"rho": 435.0}),
                         ("uptake_hill_K680", {"uptake_form": "hill"}),
                         ("rho_falls_with_glucose_full_kuhlman", {"rho_crp_coupling": "glucose"}),
                         ("rho_falls_with_glucose_weak_occ0.05", {"rho_crp_coupling": "glucose", "crp_occ_ref": 0.05}),
                         ("rho_falls_with_glucose_paper_omega10.3", {"rho_crp_coupling": "glucose",
                                                                     "loop_omega": 10.3, "loop_L0": 3.7})]:
        q, _ = fit_glucose(dict(rn.DEFAULT_PARAMS, **extra))
        q_tab = threshold_table(q)
        sens[label] = {"params": {k: q[k] for k in ["rho", "k_up", "crp_K", "excl_K", "uptake_form",
                                                    "rho_crp_coupling", "crp_occ_ref", "loop_omega"]},
                       "rho_at_1mM_glucose": float(rn.rho_eff(1000.0, q)),
                       "on_rmse_log10": float(np.sqrt(np.mean(q_tab.on_log10_err ** 2))),
                       "off_rmse_log10": float(np.sqrt(np.mean(q_tab.off_log10_err ** 2)))}
    # Inverse check: what rho does Fig 2c imply on its own?
    q_free, _ = fit_glucose(dict(rn.DEFAULT_PARAMS), fit_rho=True)

    rho_tab.to_csv(outdir / "rho_estimates_rnaseq.csv", index=False)
    tab.to_csv(outdir / "threshold_fit_vs_ozbudak.csv", index=False)
    fitted = {
        "params": {k: (v if (v is None or isinstance(v, str)) else float(v)) for k, v in p.items()},
        "derived": {"Y_max_LacY_molecules": rn.y_max(p), "basal_LacY_molecules": rn.y_max(p) / p["rho"],
                    "k_rep_max_per_min": rn.k_rep_max(p), "d_y_per_min": rn.d_y(p),
                    "rho_by_condition": {"no_glucose": rn.rho_eff(0.0, p), "glucose_1mM": rn.rho_eff(1000.0, p),
                                         "crp_knockout": rn.rho_eff(0.0, p, crp_knockout=True)},
                    "switching_drives": list(rn.switching_drives(p["rho"], p["n"]))},
        "fit": {"switch_on_rmse_log10": float(np.sqrt(np.mean(tab.on_log10_err ** 2))),
                "switch_off_rmse_log10_heldout": float(np.sqrt(np.mean(tab.off_log10_err ** 2))),
                "crp_rho_coupling": p.get("rho_crp_coupling", False)},
        "rho_free_fit_to_fig2c": float(q_free["rho"]),
        "rnaseq": {"rho": rho_tab.to_dict(orient="records"),
                   "induced_lacY_over_gat_by_group": per_group.round(4).to_dict(),
                   "crp_scaling": crp},
        "sensitivity": sens,
    }
    (outdir / "fitted_params.json").write_text(json.dumps(fitted, indent=2, default=float))
    return p, fitted, tab, rho_tab, p1k, k12


if __name__ == "__main__":
    root = Path(sys.argv[1])
    out = Path(sys.argv[2])
    p, fitted, tab, rho_tab, *_ = run_all(root, out)
    print(rho_tab.round(1).to_string())
    print(tab.round(2).to_string())
    print(json.dumps({k: fitted[k] for k in ["fit", "rho_free_fit_to_fig2c"]}, indent=1))
    print({k: round(p[k], 4) for k in ["k_up", "crp_K", "excl_K"]})
