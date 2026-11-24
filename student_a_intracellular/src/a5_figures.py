"""
A5 figures + YAML finishing step. Run after a5_validate.py.

  outputs/a5/fig_distributions.png   simulated vs real per-cell counts (histograms) for the expressed circuit genes
  outputs/a5/fig_qq.png              QQ plots: PC1 scores and pooled per-cell counts, real vs simulated
  outputs/a5/fig_phase_split.png     the exp/stat split on the 17 circuit operons: real vs simulated PC1 + fold changes
  outputs/a5/fig_lac_memory.png      lac leak / switching with the A4 bursts vs the calibrated bursts
  and fills observation_model.q_g / medium_scale into outputs/calibrated_circuit_model.yaml
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

ROOT = Path(__file__).parent.parent
A5 = ROOT / "outputs" / "a5"
REAL, SIM, OLD = "#2a78d6", "#eb6834", "#8a8984"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
                     "font.size": 10, "axes.edgecolor": "#c8c7c2", "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "bold",
                     "axes.titlesize": 10.5, "legend.frameon": False})
CRP_UP = {"lacZYA", "hflC", "tdcB", "uhpT", "cas3", "rbs", "yahN", "mhpF", "flhDC", "araJ", "nadC", "galK", "malK", "preA", "idnO", "fecE", "gyrA"}
GROWTH = {"fecE", "gyrA", "rbs", "malK"}   # fall in stationary phase for reasons outside the circuit (growth programs)


def bars(ax, real, sim, title):
    x = np.arange(5)
    r, s = np.array(real, float), np.array(sim, float)
    r, s = r / r.sum(), s / s.sum()
    ax.bar(x - 0.2, r, 0.38, color=REAL, label="PETRI-seq cells")
    ax.bar(x + 0.2, s, 0.38, color=SIM, label="simulated cells")
    ax.set_yscale("log"); ax.set_ylim(5e-5, 1.5)
    ax.set_xticks(x, ["0", "1", "2", "3", "4+"]); ax.set_title(title, loc="left")


def fig_distributions(df):
    genes = ["hflC", "rbs", "fecE", "galK", "gyrA", "lacZYA"]
    fig, axes = plt.subplots(2, 6, figsize=(16, 5.6), sharey=True)
    d = df[df.model == "calibrated"]
    for j, g in enumerate(genes):
        for i, ph in enumerate(("exp", "stat")):
            rw = d[(d.operon == g) & (d.phase == ph)].iloc[0]
            bars(axes[i, j], json.loads(rw.real_hist), json.loads(rw.sim_hist),
                 f"{g} · {'growing' if ph == 'exp' else 'stationary'}")
            axes[i, j].text(0.97, 0.95, f"mean real {rw.real_mean:.3f}\nmean sim {rw.sim_mean:.3f}", transform=axes[i, j].transAxes,
                            ha="right", va="top", fontsize=8, color=INK2)
    for i in range(2):
        axes[i, 0].set_ylabel("fraction of cells (log)")
    for j in range(6):
        axes[1, j].set_xlabel("UMIs in the cell")
    axes[0, 0].legend(loc="lower left", fontsize=8)
    fig.suptitle("Per-cell transcript counts: real PETRI-seq cells vs the calibrated simulator  (top: growing = fitted level, shape is a test; "
                 "bottom: stationary = predicted)", x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(); fig.savefig(A5 / "fig_distributions.png", dpi=170); plt.close(fig)


def hist_quantiles(counts, edges, qs):
    c = np.cumsum(counts) / np.sum(counts)
    return np.interp(qs, np.r_[0, c], edges)


def fig_qq(S, V, df):
    sc = np.load(A5 / "pc1_scores.npz")
    real = S["phase_split"]["model17_operons_pc1"]["pc1_hist_24bins"]
    lo, hi = real["range"]; edges = np.linspace(lo, hi, 25)
    qs = np.linspace(0.01, 0.99, 99)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for ph, col, nm in (("exp", REAL, "growing"), ("stat", SIM, "stationary")):
        rq = hist_quantiles(real[ph], edges, qs)
        for model, ls, mk in (("calibrated", "-", "o"),):
            s = sc[f"{model}_{ph}"]
            # orient both PC1s so that 'growing' is positive
            sq = np.quantile(s, qs)
            rqz = (rq - np.median(rq)) / (np.quantile(rq, .9) - np.quantile(rq, .1) + 1e-9)
            sqz = (sq - np.median(sq)) / (np.quantile(sq, .9) - np.quantile(sq, .1) + 1e-9)
            axes[0].plot(rqz, sqz, mk, ms=3.5, color=col, label=nm)
    lim = axes[0].get_xlim(); axes[0].plot(lim, lim, color="#c8c7c2", lw=1, zorder=0)
    axes[0].set_title("QQ: circuit-gene PC1 score (shape, robust-scaled)", loc="left")
    axes[0].set_xlabel("real quantile"); axes[0].set_ylabel("simulated quantile"); axes[0].legend()

    # QQ of pooled counts per gene: P(count >= k) real vs sim, every gene, both phases
    d = df[df.model == "calibrated"]
    for ph, col in (("exp", REAL), ("stat", SIM)):
        for _, rw in d[d.phase == ph].iterrows():
            r, s = np.array(json.loads(rw.real_hist), float), np.array(json.loads(rw.sim_hist), float)
            rs, ss = np.cumsum(r[::-1])[::-1][1:] / r.sum(), np.cumsum(s[::-1])[::-1][1:] / s.sum()
            ok = (rs > 0) & (ss > 0)
            axes[1].plot(rs[ok], ss[ok], "o", ms=3.5, color=col, alpha=0.8)
    axes[1].set_xscale("log"); axes[1].set_yscale("log")
    axes[1].plot([1e-4, 1], [1e-4, 1], color="#c8c7c2", lw=1, zorder=0)
    axes[1].set_title("Tail probabilities P(UMI ≥ k), every gene, k = 1..4", loc="left")
    axes[1].set_xlabel("real"); axes[1].set_ylabel("simulated")
    axes[1].plot([], [], "o", color=REAL, label="growing (k≥1..4)"); axes[1].plot([], [], "o", color=SIM, label="stationary")
    axes[1].legend(loc="upper left")

    # Fano: real vs sim, expressed genes, calibrated vs A3 bursts
    for model, col, nm in (("uncalibrated", OLD, "A3/A4 noise (no extrinsic)"), ("calibrated", SIM, "A5 calibrated")):
        dd = df[(df.model == model) & (df.phase == "exp") & (df.real_mean >= 0.02)]
        axes[2].plot(dd.real_fano, dd.sim_fano, "o", ms=6, color=col, label=nm)
    for _, rw in df[(df.model == "calibrated") & (df.phase == "exp") & (df.real_mean >= 0.02)].iterrows():
        axes[2].annotate(rw.operon, (rw.real_fano, rw.sim_fano), xytext=(4, 3), textcoords="offset points", fontsize=8, color=INK2)
    axes[2].plot([0.95, 1.7], [0.95, 1.7], color="#c8c7c2", lw=1, zorder=0)
    axes[2].set_title("Cell-to-cell spread (Fano) in growing cells", loc="left")
    axes[2].set_xlabel("real Fano factor"); axes[2].set_ylabel("simulated Fano factor"); axes[2].legend(loc="upper left")
    fig.tight_layout(); fig.savefig(A5 / "fig_qq.png", dpi=170); plt.close(fig)


def fig_phase_split(S, V, df):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    real = S["phase_split"]["model17_operons_pc1"]["pc1_hist_24bins"]
    lo, hi = real["range"]; edges = np.linspace(lo, hi, 25); mid = (edges[1:] + edges[:-1]) / 2
    for ph, col, nm in (("exp", REAL, "growing"), ("stat", SIM, "stationary")):
        h = np.array(real[ph], float)
        axes[0].step(mid, h / h.sum(), where="mid", color=col, lw=2, label=nm)
    axes[0].set_yscale("log"); axes[0].set_title(f"Real cells: split on 17 circuit genes  (AUC {V['real']['pca17_auc']:.2f})", loc="left")
    sc = np.load(A5 / "pc1_scores.npz")
    allv = np.r_[sc["calibrated_exp"], sc["calibrated_stat"]]
    e2 = np.linspace(np.quantile(allv, .001), np.quantile(allv, .999), 25); m2 = (e2[1:] + e2[:-1]) / 2
    for ph, col, nm in (("exp", REAL, "growing"), ("stat", SIM, "stationary")):
        h = np.histogram(sc[f"calibrated_{ph}"], e2)[0].astype(float)
        axes[1].step(m2, h / h.sum(), where="mid", color=col, lw=2, label=nm)
    axes[1].set_yscale("log")
    axes[1].set_title(f"Simulated cells: same pipeline  (AUC {V['calibrated']['pca17']['auc']:.2f})", loc="left")
    for a in axes[:2]:
        a.set_xlabel("PC1 score"); a.set_ylabel("fraction of cells (log)"); a.legend()

    d = df[df.model == "calibrated"].pivot(index="operon", columns="phase", values=["real_mean", "sim_mean"])
    cap = S["mrna_umi"]["stat"]["mean"] / S["mrna_umi"]["exp"]["mean"]
    ax = axes[2]
    for op, rw in d.iterrows():
        if rw[("real_mean", "exp")] < 0.005 or rw[("real_mean", "stat")] <= 0:
            continue
        fr = np.log2(rw[("real_mean", "stat")] / rw[("real_mean", "exp")] / cap)
        fs = np.log2(rw[("sim_mean", "stat")] / rw[("sim_mean", "exp")] / cap)
        col = OLD if op in GROWTH else SIM
        ax.plot(fr, fs, "o", ms=7, color=col)
        ax.annotate(op, (fr, fs), xytext=(4, 3), textcoords="offset points", fontsize=8, color=INK2)
    ax.axhline(0, color="#e2e1dc", lw=1); ax.axvline(0, color="#e2e1dc", lw=1)
    ax.plot([-4.5, 3.5], [-4.5, 3.5], color="#c8c7c2", lw=1)
    ax.plot([], [], "o", color=SIM, label="CRP-controlled genes"); ax.plot([], [], "o", color=OLD, label="growth genes (outside circuit)")
    ax.legend(loc="lower right", fontsize=8.5)
    ax.set_title("Stationary vs growing, per gene (log2, capture-corrected)", loc="left")
    ax.set_xlabel("real change"); ax.set_ylabel("simulated change")
    fig.tight_layout(); fig.savefig(A5 / "fig_phase_split.png", dpi=170); plt.close(fig)


def fig_lac(V):
    M = V["memory_check"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    labs = ["A4 model\n(20 proteins/burst)", "A5 calibrated\n(4.2 proteins/burst)"]
    vals = [M["uncalibrated_A4"]["uninduced"]["fraction_flipped_30h"] * 100, M["calibrated"]["uninduced"]["fraction_flipped_30h"] * 100]
    b = axes[0].bar(labs, vals, color=[OLD, SIM], width=0.55)
    for r, v in zip(b, vals):
        axes[0].text(r.get_x() + r.get_width() / 2, v + 0.4, f"{v:.1f}%", ha="center", color=INK)
    axes[0].set_ylabel("% of OFF cells that flip ON within 30 h")
    axes[0].set_title("Leak at 11 µM TMG (inside the memory window)", loc="left")
    for k, col, nm in (("uncalibrated_A4", OLD, labs[0].replace("\n", " ")), ("calibrated", SIM, labs[1].replace("\n", " "))):
        sw = M[k]["sweep_from_OFF"]
        axes[1].plot([s["TMG_uM"] for s in sw], [100 * s["fraction_ON_20h_from_OFF"] for s in sw], "-o", color=col, lw=2, ms=5, label=nm)
    axes[1].set_xscale("log"); axes[1].set_xticks([6, 8, 12, 17, 24, 34, 49], ["6", "8", "12", "17", "24", "34", "49"]); axes[1].minorticks_off()
    axes[1].set_xlabel("TMG (µM)"); axes[1].set_ylabel("% ON after 20 h (start OFF)")
    axes[1].set_title("Switching up from OFF", loc="left"); axes[1].legend(loc="upper left")
    fig.tight_layout(); fig.savefig(A5 / "fig_lac_memory.png", dpi=170); plt.close(fig)


def finish_yaml(V):
    p = ROOT / "outputs" / "calibrated_circuit_model.yaml"
    y = yaml.safe_load(p.read_text())
    c = V["calibrated"]
    y["observation_model"]["q_g"] = "same for every gene: q0 (gene-specific detection folded into medium_scale)"
    y["observation_model"]["petri_ez_rich_medium"] = {
        "description": "Expression levels in EZ Rich (PETRI-seq medium) differ from the PRECISE-1K fits. To compare with PETRI-seq "
                       "pass param_overrides={'output_scale': output_scale, 'k_tx': k_tx * lac_basal_scale}. Fitted on growing cells "
                       "only; stationary-phase levels are then predictions of the circuit.",
        "output_scale": {k: round(float(v), 4) for k, v in c["medium_scale"].items()},
        "lac_basal_scale": round(float(c["lac_basal_scale"]), 4)}
    y["noise"]["validation"] = {
        "file": "outputs/a5/validation.json",
        "growing_shape_chi2_per_dof": round(c["agreement"]["exp"]["hist_chi2"] / c["agreement"]["exp"]["hist_dof"], 2),
        "pca17_auc_sim_vs_real": [round(c["pca17"]["auc"], 3), V["real"]["pca17_auc"]],
        "lac_leak_11uM_30h": [V["memory_check"]["uncalibrated_A4"]["uninduced"]["fraction_flipped_30h"],
                              V["memory_check"]["calibrated"]["uninduced"]["fraction_flipped_30h"]]}
    p.write_text(yaml.safe_dump(y, sort_keys=False, width=130))


if __name__ == "__main__":
    S = json.loads((ROOT / "data" / "scrnaseq" / "GSE141018" / "petri_2_01_summary.json").read_text())
    V = json.loads((A5 / "validation.json").read_text())
    df = pd.read_csv(A5 / "per_gene_comparison.csv")
    for c in ("real_hist", "sim_hist"):
        df[c] = df[c].map(lambda s: json.dumps([int(float(x)) for x in s.strip("[]").replace(",", " ").split()]))
    fig_distributions(df); fig_qq(S, V, df); fig_phase_split(S, V, df); fig_lac(V); finish_yaml(V)
    print("ok")
