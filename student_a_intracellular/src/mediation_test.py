"""
A2: causal-mediation test of CRP -> target edges on PRECISE-1K (gbmMINER Fig. 5B-D logic, applied to CRP).

    perturbation of CRP  ──►  CRP activity (Crp-1 iModulon)  ──►  target gene
          └──────────── direct path (should vanish if the effect is mediated) ─────────┘

Samples: the PRECISE-1K `crp` project only (one batch): crp knockouts (fructose, glucose, glycerol; n = 9),
activating-region deletions dAr1, dAr2, dAr1dAr2 (glycerol; n = 8), wild type (fructose 1, glycerol 2).
Carbon source is a covariate in every model.

For each target Y (log2 TPM):
  total model      Y ~ D + carbon                     D = 1 for any CRP perturbation
  sign check       Y ~ KO + DD + carbon               knockout and domain deletion must move Y the same way
  conditioned      Y ~ D + C + carbon                 C = Crp-1 activity, recomputed with Y's gene left out
Causal (rule fixed before looking at results):
  BH FDR of D in the total model < 0.05, KO and DD coefficients same sign,
  p(D | C) > 0.05, p(C) < 0.05, and (beta_D - beta_D|C) / beta_D >= 0.5.

Mediator: CRP mRNA cannot be used (knockout label and crp mRNA correlate at r = -0.998; domain deletions keep normal
crp mRNA), so CRP activity is used. Leaving the tested gene out of the projection stops a gene from explaining itself.

Usage: python src/mediation_test.py   -> outputs/a2/mediation_results.csv, outputs/a2/mediation_summary.json
"""
from __future__ import annotations

import gzip
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from dima import P1K_DIR, P1K_REFS, bh  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "outputs" / "a2"
LAC = {"lacZ": "b0344", "lacY": "b0343", "lacA": "b0342", "lacI": "b0345"}


def load():
    raw = json.load(gzip.open(P1K_DIR / "p1k_precise1k.json.gz"))
    rd = lambda k: pd.read_json(io.StringIO(raw[k]))  # noqa: E731
    M, A, trn, gt = rd("M"), rd("A"), rd("trn"), rd("gene_table")
    M.attrs["thresholds"] = raw["thresholds"]
    md = pd.read_csv(P1K_DIR / "p1k_metadata_qc.csv", index_col=0)
    lt = pd.read_csv(P1K_DIR / "p1k_log_tpm_qc.csv", index_col=0)
    X = lt.sub(lt[P1K_REFS].mean(axis=1), axis=0)            # centred on the PRECISE-1K reference (= log_tpm_norm)
    return M, A, trn, gt, md, lt, X


def design(md):
    s = md[md.project == "crp"].copy()
    desc = s["Strain Description"].astype(str)
    s["KO"] = desc.str.contains("del_crp").astype(int)
    s["DD"] = desc.str.contains("crp_delAr").astype(int)
    s["D"] = s.KO | s.DD
    s["carbon"] = s["Carbon Source (g/L)"].astype(str).str.extract(r"^([a-z]+)")[0]
    return s


def ols(y, Xd):
    """OLS with t-test p-values. Xd includes the intercept column."""
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    resid = y - Xd @ beta
    dof = len(y) - Xd.shape[1]
    s2 = resid @ resid / dof
    se = np.sqrt(np.diag(s2 * np.linalg.pinv(Xd.T @ Xd)))
    t = beta / se
    return beta, 2 * stats.t.sf(np.abs(t), dof)


def crp_activity_leave_one_out(M, X, samples, gene):
    genes = M.index.intersection(X.index).drop(gene, errors="ignore")
    m_inv = np.linalg.pinv(M.loc[genes].values)
    row = list(M.columns).index("Crp-1")
    return m_inv[row] @ X.loc[genes, samples].values


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    M, A, trn, gt, md, lt, X = load()
    s = design(md)
    samples = s.index.tolist()
    carb = pd.get_dummies(s.carbon, drop_first=True).astype(float).values
    one = np.ones((len(s), 1))
    D, KO, DD = s.D.values[:, None].astype(float), s.KO.values[:, None].astype(float), s.DD.values[:, None].astype(float)

    # check the leave-nothing-out projection against the published activities
    genes = M.index.intersection(X.index)
    crp1_full = np.linalg.pinv(M.loc[genes].values)[list(M.columns).index("Crp-1")] @ X.loc[genes, samples].values
    r_full = float(np.corrcoef(crp1_full, A.loc["Crp-1", samples])[0, 1])
    # collinearity of the mRNA mediator
    crp_mrna = lt.loc["b3357", samples].values
    r_ko_mrna = float(np.corrcoef(s.KO.values, crp_mrna)[0, 1])

    reg = trn[trn.regulator.astype(str).str.lower() == "crp"]
    targets = sorted(set(reg.gene_id.astype(str)) | set(LAC.values()))
    targets = [g for g in targets if g in lt.index]
    name = gt.gene_name.to_dict()
    rows = []
    for g in targets:
        y = lt.loc[g, samples].values
        if np.allclose(y, y[0]):
            continue
        b_tot, p_tot = ols(y, np.hstack([one, D, carb]))
        b_sep, _ = ols(y, np.hstack([one, KO, DD, carb]))
        C = crp_activity_leave_one_out(M, X, samples, g)[:, None]
        b_c, p_c = ols(y, np.hstack([one, D, C, carb]))
        b_m, p_m = ols(y, np.hstack([one, D, crp_mrna[:, None], carb]))     # mRNA mediator (for the record)
        rows.append({"b_number": g, "gene": name.get(g, g), "total_effect": b_tot[1], "p_total": p_tot[1],
                     "ko_effect": b_sep[1], "dd_effect": b_sep[2], "direct_effect_given_C": b_c[1], "p_direct_given_C": p_c[1],
                     "mediator_coef": b_c[2], "p_mediator": p_c[2],
                     "frac_mediated": (b_tot[1] - b_c[1]) / b_tot[1] if b_tot[1] != 0 else np.nan,
                     "p_direct_given_mRNA": p_m[1], "p_mRNA": p_m[2],
                     "frac_mediated_mRNA": (b_tot[1] - b_m[1]) / b_tot[1] if b_tot[1] != 0 else np.nan,
                     "crp1_member": bool(g in M.index and abs(M.loc[g, "Crp-1"]) > M.attrs["thresholds"]["Crp-1"])})
    R = pd.DataFrame(rows)
    R["fdr_total"] = bh(R.p_total).values
    R["effect"] = (R.fdr_total < 0.05) & (np.sign(R.ko_effect) == np.sign(R.dd_effect))
    R["causal"] = R.effect & (R.p_direct_given_C > 0.05) & (R.p_mediator < 0.05) & (R.frac_mediated >= 0.5)
    R["causal_mRNA_mediator"] = R.effect & (R.p_direct_given_mRNA > 0.05) & (R.p_mRNA < 0.05) & (R.frac_mediated_mRNA >= 0.5)
    R = R.sort_values(["causal", "fdr_total"], ascending=[False, True])
    R.to_csv(OUT / "mediation_results.csv", index=False)
    summ = {"samples": {"knockout": int(s.KO.sum()), "domain_deletion": int(s.DD.sum()), "wild_type": int((s.D == 0).sum())},
            "targets_tested": int(len(R)), "perturbation_effect": int(R.effect.sum()), "causal": int(R.causal.sum()),
            "causal_genes": R[R.causal].gene.tolist(), "causal_that_are_Crp1_members": R[R.causal & R.crp1_member].gene.tolist(), "causal_with_mRNA_mediator": int(R.causal_mRNA_mediator.sum()),
            "r_leave_none_out_vs_published_Crp1": r_full, "r_knockout_vs_crp_mRNA": r_ko_mrna,
            "lac": R[R.gene.isin(list(LAC))][["gene", "total_effect", "fdr_total"]].round(3).to_dict("records")}
    (OUT / "mediation_summary.json").write_text(json.dumps(summ, indent=2, default=float))
    return R, summ


if __name__ == "__main__":
    R, summ = main()
    print(json.dumps(summ, indent=1, default=float))
