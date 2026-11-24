"""A6 figures (white background, plain labels). Run after a6_reproduce_aat.py and a6_dima_compare.py."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from a6_reproduce_aat import PAPER_AAT  # noqa: E402

ROOT = HERE.parent
A6 = ROOT / "outputs" / "a6"
REAL, SIM, GREY, LIGHT = "#2a78d6", "#eb6834", "#8a8984", "#d6d5d0"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 10,
                     "axes.edgecolor": "#c8c7c2", "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "bold",
                     "axes.titlesize": 10.5, "legend.frameon": False})
LABEL = {"glucose_to_glycerol": "glucose → glycerol", "glucose_to_acetate": "glucose → acetate",
         "dcrp_glucose": "Δcrp, glucose", "dcrp_glycerol": "Δcrp, glycerol", "dcyaA_LB": "ΔcyaA, LB ★",
         "dcyaA_LB_stationary": "ΔcyaA, LB stationary ★", "IPTG_50uM_glucose": "IPTG 50 µM ★",
         "IPTG_50uM_glucose_svi3": "IPTG 50 µM (2nd pair) ★"}


def fig_aat():
    r = pd.read_csv(A6 / "aat_dima_reproduction.csv", index_col=0)
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))
    a = ax[0]
    a.scatter(r.activity_1, r.activity_2, s=10, color=LIGHT, zorder=1)
    s = r[r.significant]
    a.scatter(s.activity_1, s.activity_2, s=28, color=REAL, zorder=2)
    for k, rw in s.iterrows():
        a.annotate(k, (rw.activity_1, rw.activity_2), xytext=(4, 2), textcoords="offset points", fontsize=8, color=INK2)
    lim = [min(r.activity_1.min(), r.activity_2.min()) - 2, max(r.activity_1.max(), r.activity_2.max()) + 2]
    a.plot(lim, lim, color="#c8c7c2", lw=1, zorder=0)
    a.set_xlabel("iModulon activity at aeration onset (0 min)"); a.set_ylabel("activity 10 min after aeration")
    a.set_title("Our run of Fig. 6 DiMA (AAT vs PRECISE-1K)", loc="left")
    b = ax[1]
    p = pd.Series(PAPER_AAT)
    b.scatter(p, r.loc[p.index, "difference"], color=REAL, s=34)
    for k in p.index:
        b.annotate(k, (p[k], r.loc[k, "difference"]), xytext=(4, -9), textcoords="offset points", fontsize=8, color=INK2)
    b.plot([-18, 30], [-18, 30], color="#c8c7c2", lw=1, zorder=0)
    b.set_xlabel("activity change, authors' notebook"); b.set_ylabel("activity change, our re-implementation")
    b.set_title("11 / 11 paper DiMAs reproduced (max error 0.02)", loc="left")
    fig.tight_layout(); fig.savefig(A6 / "fig_a6_aat_reproduction.png", dpi=170); plt.close(fig)


def fig_genomewide(df, S):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.1, 1]})
    a = ax[0]
    y = np.arange(len(S))
    a.barh(y + 0.2, S.real_dimas_default, 0.38, color=REAL, label="real data")
    a.barh(y - 0.2, S.sim_dimas_default, 0.38, color=SIM, label="simulated")
    for i, (r, c) in enumerate(zip(S.real_dimas_default, S.ceiling_recovered_default)):
        a.text(r + 0.8, i + 0.2, str(r), va="center", fontsize=8, color=INK2)
        a.text(0.8, i - 0.2, f"0   (best possible from circuit genes: {c})", va="center", fontsize=8, color=INK2)
    a.set_yticks(y, [LABEL[p] for p in S.perturbation]); a.invert_yaxis()
    a.set_xlabel("significant iModulons (DiMA, FDR 0.1, |Δ| > 5)")
    a.set_title("Fig. 6 readout, genome-wide projection", loc="left"); a.legend(loc="lower right")
    b = ax[1]
    for p, d in df.groupby("perturbation", sort=False):
        b.scatter(d.ceiling_difference, d.sim_difference, s=9, alpha=0.6, label=LABEL[p])
    b.plot([-5, 3], [-5, 3], color="#c8c7c2", lw=1, zorder=0)
    b.set_xlabel("activity change from the real circuit-gene changes alone")
    b.set_ylabel("activity change from the simulator")
    rr = S.set_index("perturbation").r_sim_vs_ceiling
    b.set_title(f"Simulator vs what its genes could move (median r {rr.median():.2f})", loc="left")
    b.legend(fontsize=7, loc="upper left", ncol=1)
    fig.tight_layout(); fig.savefig(A6 / "fig_a6_genomewide.png", dpi=170); plt.close(fig)


def fig_panel(P, S):
    perts = [p for p in S.perturbation if not p.startswith("IPTG")]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6))
    for a, p in zip(axes.flat, perts):
        d = P[P.perturbation == p]
        lim = max(d.panel_real_difference.abs().max(), d.panel_sim_difference.abs().max()) * 1.15 + 1
        a.axhline(0, color="#e2e1dc", lw=1); a.axvline(0, color="#e2e1dc", lw=1)
        a.plot([-lim, lim], [-lim, lim], color="#c8c7c2", lw=1)
        for _, r in d.iterrows():
            both = r.panel_real_dima and r.panel_sim_dima and np.sign(r.panel_real_difference) == np.sign(r.panel_sim_difference)
            col = SIM if both else (REAL if r.panel_real_dima else GREY)
            a.scatter(r.panel_real_difference, r.panel_sim_difference, s=46 if col != GREY else 22, color=col, zorder=3)
            if abs(r.panel_real_difference) > 3 or r.panel_sim_dima:
                a.annotate(r.imodulon, (r.panel_real_difference, r.panel_sim_difference), xytext=(4, 3),
                           textcoords="offset points", fontsize=8, color=INK2)
        s = S[S.perturbation == p].iloc[0]
        a.set_title(f"{LABEL[p]}   recovered {s.panel_recovered}/{s.panel_real_dimas}", loc="left")
        a.set_xlim(-lim, lim); a.set_ylim(-lim, lim)
        a.set_xlabel("real activity change"); a.set_ylabel("simulated activity change")
    axes.flat[0].scatter([], [], color=SIM, s=46, label="DiMA in both (same sign)")
    axes.flat[0].scatter([], [], color=REAL, s=46, label="DiMA in real data only")
    axes.flat[0].scatter([], [], color=GREY, s=22, label="not significant")
    axes.flat[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("Circuit-gene panel readout: iModulon activity changes, real vs simulated (★ = held-out data)",
                 x=0.01, ha="left", fontsize=11.5, color=INK)
    fig.tight_layout(); fig.savefig(A6 / "fig_a6_panel_dima.png", dpi=170); plt.close(fig)


def fig_genes(G):
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    a = ax[0]
    g = G[~G.is_deleted_gene & ~G.perturbation.str.startswith("IPTG")]
    for p, d in g.groupby("perturbation", sort=False):
        a.scatter(d.real_log2fc, d.sim_log2fc, s=16, alpha=0.75, label=LABEL[p])
    a.plot([-8, 4], [-8, 4], color="#c8c7c2", lw=1); a.axhline(0, color="#e2e1dc", lw=1); a.axvline(0, color="#e2e1dc", lw=1)
    sl = np.polyfit(g.real_log2fc, g.sim_log2fc, 1)[0]
    a.set_title(f"Circuit genes: right direction, ~{1/sl:.1f}× too small", loc="left")
    a.set_xlabel("real log2 fold change"); a.set_ylabel("simulated log2 fold change"); a.legend(fontsize=7)
    for gname in ("tdcB", "malK", "preA"):
        r = g[(g.gene == gname) & (g.perturbation == "dcyaA_LB")]
        if len(r):
            a.annotate(f"{gname} (ΔcyaA)", (r.real_log2fc.iloc[0], r.sim_log2fc.iloc[0]), xytext=(5, 4),
                       textcoords="offset points", fontsize=8, color=INK2)
    b = ax[1]
    L = G[G.perturbation.str.startswith("IPTG") & G.gene.isin(["lacZ", "lacY", "lacA"])]
    x = np.arange(3)
    for i, (p, d) in enumerate(L.groupby("perturbation", sort=False)):
        d = d.set_index("gene").loc[["lacZ", "lacY", "lacA"]]
        b.bar(x - 0.27 + i * 0.27, d.real_log2fc, 0.25, yerr=d.real_log2fc_sd, color=REAL if i == 0 else "#7fb0ea",
              label=f"real, {'pair 1' if i == 0 else 'pair 2'}", capsize=2)
    b.bar(x + 0.27, L.groupby("gene").sim_log2fc.mean().loc[["lacZ", "lacY", "lacA"]], 0.25, color=SIM, label="simulated")
    b.set_xticks(x, ["lacZ", "lacY", "lacA"]); b.set_ylabel("log2 fold change, +50 µM IPTG")
    gap = L.groupby("perturbation").apply(lambda d: d.sim_log2fc.mean() - d.real_log2fc.mean())
    b.set_title(f"Lac induction: model {2 ** gap.min():.0f}–{2 ** gap.max():.0f}× stronger than these strains", loc="left")
    b.legend()
    fig.tight_layout(); fig.savefig(A6 / "fig_a6_gene_level.png", dpi=170); plt.close(fig)


def fig_hist(df, P):
    """Fig. 6 'activity histogram': P1K distribution with the matched real and simulated samples (Maltose, dcyaA_LB)."""
    sys.path.insert(0, str(HERE))
    from dima import KnowledgeBase
    kb = KnowledgeBase.load(with_log_tpm=False)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for a, (im, p) in zip(ax, [("Maltose", "dcyaA_LB"), ("Crp-1", "dcrp_glycerol")]):
        a.hist(kb.A.loc[im], bins=60, color=LIGHT, label="all PRECISE-1K samples")
        r = df[(df.perturbation == p) & (df.imodulon == im)].iloc[0]
        q = P[(P.perturbation == p) & (P.imodulon == im)].iloc[0]
        ymax = a.get_ylim()[1]
        for x, c, lab, ls in [(r.real_activity_control, GREY, "real control", "-"), (r.real_activity_perturbed, REAL, "real perturbed", "-"),
                              (r.sim_activity_perturbed, SIM, "simulated, genome-wide readout", "--"),
                              (r.real_activity_control + q.panel_sim_difference, SIM, "simulated, circuit-gene panel", "-")]:
            a.axvline(x, color=c, lw=2, ls=ls, label=lab)
        a.set_title(f"{im} activity · {LABEL[p]}", loc="left"); a.set_xlabel(f"{im} iModulon activity"); a.set_ylabel("samples")
        a.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout(); fig.savefig(A6 / "fig_a6_activity_histograms.png", dpi=170); plt.close(fig)


if __name__ == "__main__":
    df = pd.read_csv(A6 / "perturbation_dima_comparison.csv")
    S = pd.read_csv(A6 / "perturbation_summary.csv")
    P = pd.read_csv(A6 / "panel_dima_comparison.csv")
    G = pd.read_csv(A6 / "gene_level_comparison.csv")
    fig_aat(); fig_genomewide(df, S); fig_panel(P, S); fig_genes(G); fig_hist(df, P)
    print("ok")
