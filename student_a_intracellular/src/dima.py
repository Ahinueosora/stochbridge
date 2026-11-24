"""
PRECISE-1K knowledge-base workflow (Lamoureux et al. 2023, Fig. 6 / notebook workflow__add_new_data.ipynb), re-implemented
so it runs on current pandas (pymodulon 0.2.1's loader and Index '&' break on pandas 2.x). Same maths:

  1. new log2(TPM+1) data, restricted to genes in PRECISE-1K, centred on the P1K reference (p1k_00001/p1k_00002, wt glucose M9)
  2. inferred activities  A_new = pinv(M[shared genes]) @ X_new[shared genes]
  3. DiMA (pymodulon.util.dima): for every iModulon, null = log-normal fitted to |activity difference| between replicates
     (one pair per project:condition, as pymodulon does), p = 1 - CDF(|a2 - a1|), Benjamini-Hochberg FDR, |difference| > threshold.
     Paper settings for AAT: fdr = 0.05, threshold = 9.
"""
from __future__ import annotations

import gzip
import io
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).parent.parent
P1K_DIR = ROOT / "data" / "precise1k"
P1K_REFS = ["p1k_00001", "p1k_00002"]


class KnowledgeBase:
    """M, A, tables of PRECISE-1K plus any samples projected onto it."""

    def __init__(self, M, A, sample_table, imodulon_table, gene_table, log_tpm=None):
        self.M, self.A, self.sample_table = M, A, sample_table
        self.imodulon_table, self.gene_table, self.log_tpm = imodulon_table, gene_table, log_tpm
        self._null = None

    # ------------------------------------------------------------------ loading
    @classmethod
    def load(cls, d: Path = P1K_DIR, with_log_tpm=True):
        raw = json.load(gzip.open(d / "p1k_precise1k.json.gz"))
        rd = lambda k: pd.read_json(io.StringIO(raw[k]))  # noqa: E731
        M, A = rd("M"), rd("A")
        st, it, gt = rd("sample_table"), rd("imodulon_table"), rd("gene_table")
        M.index = M.index.astype(str); A.columns = A.columns.astype(str)
        lt = pd.read_csv(d / "p1k_log_tpm_qc.csv", index_col=0) if with_log_tpm else None
        return cls(M, A, st, it, gt, lt)

    # ------------------------------------------------------------------ projection (Fig. 6 step)
    def center(self, log_tpm: pd.DataFrame):
        genes = self.log_tpm.index.intersection(log_tpm.index)
        ref = self.log_tpm.loc[genes, P1K_REFS].mean(axis=1)
        return log_tpm.loc[genes].sub(ref, axis=0)

    def infer_activities(self, centred: pd.DataFrame):
        shared = self.M.index.intersection(centred.index)
        m_inv = np.linalg.pinv(self.M.loc[shared].values)
        return pd.DataFrame(m_inv @ centred.loc[shared].values, index=self.M.columns, columns=centred.columns)

    def add_samples(self, log_tpm: pd.DataFrame, metadata: pd.DataFrame, centred=False):
        """Project new samples and append them (activities + metadata) to the knowledge base."""
        X = log_tpm if centred else self.center(log_tpm)
        A_new = self.infer_activities(X)
        self.A = pd.concat([self.A.drop(columns=A_new.columns, errors="ignore"), A_new], axis=1)
        md = metadata.loc[A_new.columns]
        self.sample_table = pd.concat([self.sample_table.drop(index=md.index, errors="ignore"), md])
        self._null = None
        return A_new

    # ------------------------------------------------------------------ DiMA
    def samples(self, s):
        if isinstance(s, str):
            proj, cond = s.split(":", 1)
            st = self.sample_table
            idx = st[(st.project == proj) & (st.condition == cond)].index.tolist()
            if not idx:
                raise ValueError(f"no samples for {s}")
            return idx
        return list(s)

    def null_cdfs(self, A=None):
        """pymodulon.util.dima null: |difference| of one replicate pair per project:condition, log-normal per iModulon."""
        if self._null is not None and A is None:
            return self._null
        A = self.A if A is None else A
        diffs = {}
        st = self.sample_table.loc[self.sample_table.index.intersection(A.columns)]
        for name, g in st.groupby(["project", "condition"]):
            for i1, i2 in combinations(g.index, 2):
                diffs[":".join(map(str, name))] = (A[i1] - A[i2]).abs()     # last pair wins, as in pymodulon
        D = pd.DataFrame(diffs)
        null = {k: stats.lognorm(*stats.lognorm.fit(D.loc[k].values)) for k in A.index}
        if A is self.A:
            self._null = null
        return null

    def dima(self, s1, s2, threshold=5.0, fdr=0.1, full=True):
        """Full table: activity in s1 and s2, difference (s2 - s1), p, BH q, significant flag."""
        l1, l2 = self.samples(s1), self.samples(s2)
        null = self.null_cdfs()
        a1, a2 = self.A[l1].mean(axis=1), self.A[l2].mean(axis=1)
        diff = a2 - a1
        p = pd.Series({k: 1 - null[k].cdf(abs(diff[k])) for k in self.A.index})
        q = bh(p)
        res = pd.DataFrame({"activity_1": a1, "activity_2": a2, "difference": diff, "pvalue": p, "qvalue": q})
        res["significant"] = (res.qvalue < fdr) & (res.difference.abs() > threshold)
        res = res.sort_values("difference", ascending=False)
        return res if full else res[res.significant]


def bh(p: pd.Series):
    """Benjamini-Hochberg q-values (pymodulon.util.FDR keeps rows with q < fdr; same ordering rule)."""
    n = len(p)
    order = p.sort_values().index
    ranked = p.loc[order].values * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    return pd.Series(np.minimum(q, 1), index=order).loc[p.index]


# ---------------------------------------------------------------------- counts -> log TPM (for K-12 modulome samples)
def log_tpm_from_counts(counts: pd.DataFrame, gene_lengths: pd.Series):
    genes = counts.index.intersection(gene_lengths.index)
    c = counts.loc[genes]
    rpk = c.div(gene_lengths.loc[genes].values / 1000.0, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    return np.log2(tpm + 1)


# ---------------------------------------------------------------------- circuit-gene panel readout (A6 extension)
class Panel:
    """iModulon activities estimated from a small gene panel (the genes a circuit model produces).

    For iModulon k with panel members C_k (genes of the panel inside the iModulon, |M_gk| > threshold):
        a_k = sum_{g in C_k} M_gk x_g / sum_{g in C_k} M_gk^2     (least squares on the member genes only)
    x is log2(TPM+1) centred on the PRECISE-1K reference. 'fidelity' = Pearson r between a_k and the full ICA activity
    over all PRECISE-1K samples: how well the panel can read that iModulon at all."""

    def __init__(self, kb: KnowledgeBase, genes, thresholds: dict, min_fidelity=0.5):
        genes = [g for g in genes if g in kb.M.index]
        th = pd.Series(thresholds)
        self.W = {}
        for im in kb.M.columns:
            mem = [g for g in genes if abs(kb.M.loc[g, im]) > th[im]]
            if mem:
                w = kb.M.loc[mem, im]
                self.W[im] = w / float((w ** 2).sum())
        Xp = kb.center(kb.log_tpm)
        Ap = self.activities(Xp)
        self.fidelity = pd.Series({im: float(np.corrcoef(Ap.loc[im], kb.A.loc[im, Ap.columns])[0, 1]) for im in Ap.index})
        self.members = {im: list(w.index) for im, w in self.W.items()}
        self.readable = [im for im in self.W if self.fidelity[im] >= min_fidelity]
        # null: |difference| of one replicate pair per project:condition in PRECISE-1K, log-normal (as pymodulon)
        st = kb.sample_table.loc[Ap.columns]
        diffs = {}
        for name, g in st.groupby(["project", "condition"]):
            for i1, i2 in combinations(g.index, 2):
                diffs[":".join(map(str, name))] = (Ap[i1] - Ap[i2]).abs()
        D = pd.DataFrame(diffs)
        self.null = {im: stats.lognorm(*stats.lognorm.fit(D.loc[im].values[D.loc[im].values > 0])) for im in self.W}

    def activities(self, centred: pd.DataFrame):
        return pd.DataFrame({im: w @ centred.loc[w.index] for im, w in self.W.items()}).T

    def dima(self, a1: pd.DataFrame, a2: pd.DataFrame, ims=None):
        ims = ims or self.readable
        d = (a2.mean(axis=1) - a1.mean(axis=1)).loc[ims]
        p = pd.Series({k: 1 - self.null[k].cdf(abs(d[k])) for k in ims})
        return d, p, bh(p)
