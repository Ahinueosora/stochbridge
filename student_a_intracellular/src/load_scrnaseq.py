"""
A5 - load_scrnaseq: PETRI-seq (Blattman et al. 2020, GEO GSE141018) count matrix -> standard cell x operon frame.

    cells, X, operons = load_petri("data/scrnaseq/GSE141018/GSM4489500_..._gene_matrix.txt.gz")
    S = summarise(cells, X, operons)      # the summary used by fit_noise_params.py

GEO file format: tab-separated, first row "Cell Barcode" + operon names ("U00096:<operon>" for E. coli,
"CP000255...:" for S. aureus, "anti-" prefix for antisense), one row per barcode ("SB442_bc1_10_bc2_11_bc3_22").
The 2.01 matrix decompresses to ~1.3 GB of text, so it is streamed line by line into a sparse matrix.

Standard frame:
    cells    pandas DataFrame, one row per kept E. coli cell: barcode, bc1, phase ("exp"/"stat"), ecoli_umi, mrna_umi
    X        scipy.sparse CSR (cells x operons), raw mRNA UMI counts
    operons  list of operon names (U00096 sense mRNA operons)

Filters (paper Methods + GEO record): E. coli cell = >90% of UMIs on U00096, <20 S. aureus UMIs, >=15 UMIs;
mRNA = U00096 sense operons minus rRNA, tRNA-only operons, ssrA, ffs, rnpB; keep cells with >=15 mRNA UMIs.
Phase (experiment 2.01): exponential cells in RT-plate columns 1-6, stationary in 7-12 -> stationary if (bc1-1) % 12 >= 6.

The summary in data/scrnaseq/GSE141018/petri_2_01_summary.json was produced by the same steps run inside the
desktop browser (the analysis workspace cannot reach NCBI); running summarise() on the downloaded file regenerates it.
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

TRNA = re.compile(r"^(ala|arg|asn|asp|cys|gln|glt|gly|his|ile|leu|lys|met|phe|pro|ser|thr|trp|tyr|val|sel)[A-Z]$")
# the model's transcribed species -> PETRI-seq operon (operon names compress gene letters, e.g. lacZYA, galETKM)
MODEL_OPERONS = {"lacZYA": r"^lacZYA$", "hflC": r"hflXKC", "tdcB": r"^tdcABCDEFG", "uhpT": r"^uhpT$", "cas3": r"^cas3$",
                 "rbs": r"^rbsD", "yahN": r"^yahN$", "fecE": r"^fecABCDE$", "mhpF": r"^mhpABCDFE$", "flhDC": r"^flhDC$",
                 "araJ": r"^araJ$", "nadC": r"^nadC$", "galK": r"^galETKM$", "gyrA": r"^gyrA$",
                 "malK": r"^malK-lamB-malM$", "preA": r"^preTA$", "idnO": r"^idnDOTR$"}
EXTRA_OPERONS = {"lacI": r"^mhpR-lacI$", "crp": r"^crp$", "cyaA": r"^cyaA$", "gat": r"^gatYZABCD$"}
# which A2 CRP output(s) each operon stands for in the model (rbs operon = rbsA/B/C/R; use rbsA; flhDC -> flhD; malK operon -> malK)
OPERON_TO_MODEL_OUTPUT = {"hflC": "hflC", "tdcB": "tdcB", "uhpT": "uhpT", "cas3": "cas3", "rbs": "rbsA", "yahN": "yahN",
                          "fecE": "fecE", "mhpF": "mhpF", "flhDC": "flhD", "araJ": "araJ", "nadC": "nadC", "galK": "galK",
                          "gyrA": "gyrA", "malK": "malK", "preA": "preA", "idnO": "idnO"}


def is_rna_operon(name: str) -> bool:
    if name == "rRNA" or re.search(r"ssrA|ffs|rnpB", name):
        return True
    return all(TRNA.match(x) for x in name.split("-"))


def phase_of_bc1(bc1: int) -> str:
    return "stat" if (bc1 - 1) % 12 >= 6 else "exp"


def load_petri(path, min_umi=15):
    """Stream the GEO matrix; return (cells DataFrame, CSR mRNA count matrix, operon names)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        is_ec = np.array([h.startswith("U00096:") for h in header])
        names = [h.split(":", 1)[1] if ":" in h else h for h in header]
        mrna = np.array([is_ec[k] and not names[k].startswith("anti-") and not is_rna_operon(names[k])
                         for k in range(len(header))])
        mcols = np.where(mrna)[0]
        col_of = {k: j for j, k in enumerate(mcols)}
        rows, cols, vals, meta = [], [], [], []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            v = np.array(parts[1:], dtype=float)
            nz = np.nonzero(v)[0]
            ec = v[nz][is_ec[nz]].sum()
            sa = v[nz][~is_ec[nz]].sum()
            tot = ec + sa
            if tot < min_umi or ec / tot <= 0.9 or sa >= 20:
                continue
            mk = nz[mrna[nz]]
            m_tot = v[mk].sum()
            if m_tot < min_umi:
                continue
            i = len(meta)
            bc1 = int(re.search(r"bc1_(\d+)_", parts[0]).group(1))
            meta.append({"barcode": parts[0], "bc1": bc1, "phase": phase_of_bc1(bc1), "ecoli_umi": ec, "mrna_umi": m_tot})
            rows += [i] * len(mk); cols += [col_of[k] for k in mk]; vals += list(v[mk])
    cells = pd.DataFrame(meta)
    X = sparse.csr_matrix((vals, (rows, cols)), shape=(len(cells), len(mcols)))
    return cells, X, [names[k] for k in mcols]


def operon_counts(X, operons, pattern):
    idx = [j for j, n in enumerate(operons) if re.search(pattern, n)]
    return np.asarray(X[:, idx].sum(1)).ravel() if idx else np.zeros(X.shape[0])


def auc_and_balacc(score, is_exp):
    from sklearn.metrics import roc_auc_score
    auc = roc_auc_score(is_exp, score)
    auc = max(auc, 1 - auc)
    o = np.argsort(score); e = is_exp[o]; nE, nS = e.sum(), (~e).sum()
    cE = np.cumsum(e); cS = np.cumsum(~e)
    bal = np.maximum((cS / nS + (nE - cE) / nE) / 2, (cE / nE + (nS - cS) / nS) / 2).max()
    return float(auc), float(bal)


def paper_pca(X, mrna_umi, min_operon_counts=5, n_components=2):
    """Blattman et al. PCA: drop operons with <= 5 counts, normalise to cell total x geometric mean, log(x+1),
    z-score per operon, PCA. Returns (scores, loadings, kept operon indices)."""
    from sklearn.decomposition import PCA
    keep = np.where(np.asarray(X.sum(0)).ravel() > min_operon_counts)[0]
    Xk = X[:, keep].toarray().astype(float)
    gm = np.exp(np.mean(np.log(mrna_umi)))
    Z = np.log(Xk / mrna_umi[:, None] * gm + 1)
    Z = (Z - Z.mean(0)) / np.where(Z.std(0) > 0, Z.std(0), 1)
    p = PCA(n_components=n_components, random_state=0).fit(Z)
    return p.transform(Z), p.components_, keep


def counts_pca(C, mrna_umi, gm=None):
    """Same normalisation on a small cells x genes count matrix (C), 1 component; sign: fecE column positive if present."""
    from sklearn.decomposition import PCA
    gm = gm or np.exp(np.mean(np.log(mrna_umi)))
    Z = np.log(C / mrna_umi[:, None] * gm + 1)
    sd = Z.std(0)
    Z = (Z - Z.mean(0)) / np.where(sd > 0, sd, 1)
    p = PCA(n_components=1, random_state=0).fit(Z)
    return p.transform(Z)[:, 0], p.components_[0]


def summarise(cells, X, operons):
    """Rebuild data/scrnaseq/GSE141018/petri_2_01_summary.json from the parsed matrix."""
    is_exp = (cells.phase == "exp").values
    tot = cells.mrna_umi.values.astype(float)
    S = {"cells": {"n_cells_pca": int(len(cells)), "n_exp": int(is_exp.sum()), "n_stat": int((~is_exp).sum())}}
    sc, load, keep = paper_pca(X, tot)
    S["cells"]["n_operons_pca"] = int(len(keep))
    S["phase_split"] = {"all_operons_pc1": dict(zip(["auc", "balanced_accuracy"], auc_and_balacc(sc[:, 0], is_exp)))}
    S["phase_split"]["total_mrna_umi_only"] = dict(zip(["auc", "balanced_accuracy"], auc_and_balacc(tot, is_exp)))
    S["mrna_umi"] = {ph: {"mean": float(tot[m].mean()),
                          "quantiles": {str(q): float(np.quantile(tot[m], q)) for q in (.05, .1, .25, .5, .75, .9, .95)}}
                     for ph, m in (("exp", is_exp), ("stat", ~is_exp))}
    gm = np.exp(np.mean(np.log(tot)))
    ops = {}
    for g, pat in {**MODEL_OPERONS, **EXTRA_OPERONS}.items():
        v = operon_counts(X, operons, pat)
        rec = {}
        for ph, m in (("exp", is_exp), ("stat", ~is_exp)):
            x = v[m]
            rec[ph] = [float(x.mean()), float(x.var(ddof=1)), [int((x == k).sum()) for k in range(4)] + [int((x >= 4).sum())],
                       float((x / tot[m] * gm).mean())]
        ops[g] = rec
    S["operons"] = ops
    C = np.column_stack([operon_counts(X, operons, p) for p in MODEL_OPERONS.values()])
    s1, w = counts_pca(C, tot, gm)
    S["phase_split"]["model17_operons_pc1"] = dict(zip(["auc", "balanced_accuracy"], auc_and_balacc(s1, is_exp)))
    S["phase_split"]["model17_operons_pc1"]["loadings"] = dict(zip(MODEL_OPERONS, map(float, w)))
    return S


if __name__ == "__main__":
    cells, X, operons = load_petri(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/scrnaseq/GSE141018/petri_2_01_summary_regenerated.json")
    out.write_text(json.dumps(summarise(cells, X, operons), indent=1))
    print(f"{len(cells)} cells x {len(operons)} operons -> {out}")
