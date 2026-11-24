"""
A4: reproduce Vilar, Guet & Leibler 2003 (J Cell Biol 161:471) Fig 3 with the A3 v3 circuit, plus the SSA hysteresis
check. Vilar's parameters are NOT used; the point is that our SSA gives all-or-none (bimodal) induction.

    simulate()  -> outputs/a4/*.parquet, outputs/a4/vilar_results.json, outputs/summary_stats.json
    figures()   -> outputs/a4/fig_a4_*.png   (also called from notebooks/03_vilar_reproduction.ipynb)

Vilar Fig 3 panels and our analogues (doses scaled to our model's memory zone 4.1-30 uM):
    3a  sub-optimal inducer: single cells switch abruptly at random times; 2000-cell average is smooth;
        population curve (induced cells grow slower) lies below the average      -> T_sub (~50% ON by 20 h)
    3b  maintenance dose inside the memory zone, from both histories             -> T_maint = 11 uM
    3c  high inducer: everyone switches, average ~ population                    -> 500 uM
Molecular level = promoter flips and mRNA bursts (exact SSA); cellular = one cell's LacY; population = average and
growth-weighted population curve.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import reaction_network as rn  # noqa: E402
from io_schema import Trajectory, bimodality, bimodality_of_counts  # noqa: E402
from perturbation import InducerLevelChange, OperatorDecoy  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "a4"
YAML = ROOT / "outputs" / "circuit_model.yaml"
GROWTH_RATIO_INDUCED = 0.85     # assumption: induced cells grow 15% slower (lac burden; Vilar used measured rates)
T_MAINT = 11.0                  # geometric middle of the 4.1-30 uM memory zone
T_HIGH = 500.0


def first_switch_times(tr: Trajectory):
    """Time each cell first crosses the ON/OFF threshold (inf if never)."""
    W = tr.wide("Y"); thr = tr.meta["on_threshold_Y"]; t = W.columns.values
    up = tr.meta["perturbation"]["initial_state"] == "uninduced"
    hit = (W.values > thr) if up else (W.values < thr)
    first = np.where(hit.any(1), t[np.argmax(hit, 1)], np.inf)
    return first


def population_curve(tr: Trajectory, ratio=GROWTH_RATIO_INDUCED):
    """Growth-weighted population mean LacY: a cell's lineage grows at lambda before it switches ON and
    ratio * lambda after (Vilar's population level). Returns times, cell average, population average."""
    lam = np.log(2) / 90.0
    W = tr.wide("Y"); t = W.columns.values
    ts = first_switch_times(tr)
    tt = t[None, :]
    on_time = np.clip(tt - ts[:, None], 0, None)
    logw = lam * tt - (1 - ratio) * lam * on_time
    w = np.exp(logw - logw.max(0, keepdims=True))
    return t, W.values.mean(0), (w * W.values).sum(0) / w.sum(0)


# ======================================================================= simulate
def simulate(n=1000, n_traj=2000):
    from ssa_simulator import StochasticCircuit
    OUT.mkdir(parents=True, exist_ok=True)
    sc = StochasticCircuit(YAML, mode="fast", include_outputs=True)
    p = sc.params
    R = {"engine": {"mode": "fast", "compile_s": sc.compile_s}, "growth_ratio_induced_assumed": GROWTH_RATIO_INDUCED}

    # ---- 1. SSA hysteresis: fraction ON after 20 h, 1000 cells per history per dose
    Ts = np.round(np.geomspace(1, 200, 16), 2)
    rows = []
    for T in Ts:
        for hist in ("uninduced", "induced"):
            tr = sc.run(InducerLevelChange(float(T), initial_state=hist, n_cells=n, seed=int(T * 100) + (hist == "induced")))
            snap = tr.at(tr.times[-1]).Y_count.values
            rows.append({"TMG_uM": float(T), "history": hist, "fraction_ON_20h": float((snap > tr.meta["on_threshold_Y"]).mean()),
                         "mean_Y": float(snap.mean()), "median_Y": float(np.median(snap)),
                         "dip_p": bimodality_of_counts(snap)["dip_p"], "wall_s": tr.meta["wall_time_s"]})
    hyst = pd.DataFrame(rows); hyst.to_csv(OUT / "ssa_hysteresis.csv", index=False)

    def crossing(hist, level):
        d = hyst[hyst.history == hist].sort_values("TMG_uM")
        x, f = np.log10(d.TMG_uM.values), d.fraction_ON_20h.values
        if hist == "induced":       # switching OFF happens as T decreases: fraction still ON falls below 1 - level
            f = 1 - f
            x, f = x[::-1], f[::-1]
        i = np.where(f >= level)[0]
        if not len(i) or i[0] == 0:
            return None
        j = i[0]
        return float(10 ** np.interp(level, [f[j - 1], f[j]], [x[j - 1], x[j]]))
    an_off, an_on = rn.thresholds_tmg(0.0, p)
    R["ssa_hysteresis"] = {
        "switch_on_T50_uM": crossing("uninduced", 0.5), "switch_on_T95_uM": crossing("uninduced", 0.95),
        "switch_off_T50_uM": crossing("induced", 0.5), "switch_off_T95_uM": crossing("induced", 0.95),
        "analytic_switch_on_uM": an_on, "analytic_switch_off_uM": an_off,
        "ode_20h_switch_on_uM": 36.4, "ode_20h_switch_off_uM": 2.94, "ozbudak_paper_uM": [3, 30], "n_cells_per_point": n}

    # ---- 2. Vilar 3a: sub-optimal dose where ~half the OFF cells switch within 20 h
    T_sub = R["ssa_hysteresis"]["switch_on_T50_uM"] or 40.0
    T_sub = float(np.round(T_sub, 1))
    R["T_sub_uM"] = T_sub
    runs = {"vilar3a_sub": InducerLevelChange(T_sub, n_cells=n_traj, t_end_min=1800, seed=7, name="vilar3a_sub"),
            "vilar3c_high": InducerLevelChange(T_HIGH, n_cells=n_traj, t_end_min=1800, seed=8, name="vilar3c_high"),
            "vilar3b_maint_from_OFF": InducerLevelChange(T_MAINT, n_cells=n, t_end_min=1800, seed=9, name="vilar3b_maint_from_OFF"),
            "vilar3b_maint_from_ON": InducerLevelChange(T_MAINT, initial_state="induced", n_cells=n, t_end_min=1800, seed=10,
                                                        name="vilar3b_maint_from_ON")}
    R["vilar"] = {}
    for key, pt in runs.items():
        tr = sc.run(pt)
        tr.to_parquet(OUT / f"trajectory_{key}.parquet")
        s = tr.write_summary(OUT / f"summary_stats_{key}.json")
        ts = first_switch_times(tr)
        t, avg, pop = population_curve(tr)
        sw = ts[np.isfinite(ts)]
        R["vilar"][key] = {"TMG_uM": pt.TMG_uM, "history": pt.initial_state, "n_cells": pt.n_cells,
                           "fraction_switched_30h": float(np.isfinite(ts).mean()),
                           "switch_time_median_min": float(np.median(sw)) if len(sw) else None,
                           "switch_time_iqr_min": [float(np.percentile(sw, 25)), float(np.percentile(sw, 75))] if len(sw) > 3 else None,
                           "final_cell_avg_Y": float(avg[-1]), "final_population_avg_Y": float(pop[-1]),
                           "summary_stats": s}
        if key == "vilar3a_sub":            # the canonical Module 1 example output (brief: outputs/summary_stats.json)
            (ROOT / "outputs" / "summary_stats.json").write_text(json.dumps(s, indent=2, default=float))
            # how fast is one switch? time from 10% to 90% of the ON level, per switching cell
            W = tr.wide("Y").values; tg = tr.times; yon = tr.meta["Y_on"]
            dur = []
            for row in W:
                a = np.where(row > 0.1 * yon)[0]; b = np.where(row > 0.9 * yon)[0]
                if len(a) and len(b) and b[0] > a[0]:
                    dur.append(tg[b[0]] - tg[a[0]])
            R["vilar"][key]["single_cell_rise_10_90_median_min"] = float(np.median(dur)) if dur else None
            R["vilar"][key]["population_rise_10_90_min"] = float(
                np.interp(0.9 * avg[-1], avg, tg) - np.interp(0.1 * avg[-1], avg, tg))

    # ---- 3. bimodal vs noisy-unimodal: wild type at T_sub vs graded rho = 5 strain at a dose with the same mean
    tr = Trajectory.from_parquet(OUT / "trajectory_vilar3a_sub.parquet")
    hist_times = [0, 240, 480, 720, 1200, 1800]
    R["bimodality_wt"] = {int(t): bimodality_of_counts(tr.at(t).Y_count.values) | {"fraction_ON": float((tr.at(t).Y_count.values > tr.meta["on_threshold_Y"]).mean())}
                          for t in hist_times}
    q5 = dict(p, rho=5.0)
    Tg = np.geomspace(0.1, 300, 200)
    y5 = np.array([rn.simulate(T, 0.0, q5, induced=False) for T in Tg])
    T5 = float(np.interp(np.log10(tr.at(1800).Y_count.mean()), np.log10(y5), Tg))
    tr5 = sc.run(OperatorDecoy(5, TMG_uM=T5, n_cells=n_traj, t_end_min=1800, seed=11, name="rho5_graded"))
    tr5.to_parquet(OUT / "trajectory_rho5_graded.parquet")
    R["rho5_control"] = {"TMG_uM": T5, "matched_mean_Y": float(tr.at(1800).Y_count.mean()),
                         "rho5_mean_Y": float(tr5.at(1800).Y_count.mean()),
                         "bimodality_rho5": {int(t): bimodality_of_counts(tr5.at(t).Y_count.values) for t in hist_times}}

    # ---- 4. k_unrep sensitivity (noise knob; rho fixed so the mean model is unchanged)
    R["k_unrep_sensitivity"] = []
    for ku in (0.02, 0.1, 0.5):
        for T, hist in ((T_MAINT, "uninduced"), (T_sub, "uninduced"), (T_MAINT, "induced")):
            trk = sc.run(InducerLevelChange(T, initial_state=hist, n_cells=n, seed=21, param_overrides={"k_unrep": ku}))
            snap = trk.at(trk.times[-1]).Y_count.values
            R["k_unrep_sensitivity"].append({"k_unrep": ku, "TMG_uM": T, "history": hist,
                                             "fraction_ON_20h": float((snap > trk.meta["on_threshold_Y"]).mean())})
    (OUT / "vilar_results.json").write_text(json.dumps(R, indent=2, default=float))
    return R


# ======================================================================= figures
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, RED, GREEN, PURPLE = "#2a78d6", "#eb6834", "#c4314b", "#2e7d32", "#7b3fa0"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
                         "figure.facecolor": "white", "axes.facecolor": "white"})
    return plt


def fig_vilar3(plt=None, save=True):
    plt = plt or _style()
    R = json.loads((OUT / "vilar_results.json").read_text())
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6))
    for ax, key, title in [(axs[0], "vilar3a_sub", f"a. Sub-optimal TMG ({R['T_sub_uM']:g} µM)"),
                           (axs[2], "vilar3c_high", f"c. High TMG ({T_HIGH:g} µM)")]:
        tr = Trajectory.from_parquet(OUT / f"trajectory_{key}.parquet")
        ts = first_switch_times(tr)
        W = tr.wide("Y"); h = W.columns.values / 60
        order = np.argsort(ts)
        finite = order[np.isfinite(ts[order])]
        picks = finite[[int(len(finite) * q) for q in (0.15, 0.5, 0.85)]] if len(finite) > 10 else finite[:3]
        for c, col in zip(picks, [ORANGE, GREEN, BLUE]):
            ax.plot(h, W.values[c], color=col, lw=1.2, alpha=0.9)
        t, avg, pop = population_curve(tr)
        ax.plot(t / 60, avg, color=INK, lw=2.6, label=f"average of {tr.n_cells:,} cells")
        ax.plot(t / 60, pop, color=INK, lw=1.8, ls="--", label="population (ON cells grow 15% slower)")
        ax.set_title(title, loc="left", fontsize=12); ax.set_xlabel("hours after adding TMG")
        ax.set_ylabel("LacY per cell")
        ax.legend(loc="upper left" if key == "vilar3a_sub" else "lower right", fontsize=9)
        if key == "vilar3a_sub":
            ax.set_ylim(top=ax.get_ylim()[1] * 1.18)
    ax = axs[1]
    for key, col, lab in [("vilar3b_maint_from_ON", ORANGE, "cells that start ON"),
                          ("vilar3b_maint_from_OFF", BLUE, "cells that start OFF")]:
        tr = Trajectory.from_parquet(OUT / f"trajectory_{key}.parquet")
        W = tr.wide("Y"); h = W.columns.values / 60
        for row in W.values[:40]:
            ax.semilogy(h, row + 1, color=col, lw=0.6, alpha=0.35)
        ax.semilogy(h, W.values.mean(0) + 1, color=col, lw=2.6, label=f"{lab}: {100 * R['vilar'][key]['fraction_switched_30h']:.1f}% flip in 30 h")
    ax.set_title(f"b. Maintenance dose inside the memory zone ({T_MAINT:g} µM)", loc="left", fontsize=12)
    ax.set_xlabel("hours"); ax.set_ylabel("LacY per cell + 1 (log)")
    ax.set_ylim(0.8, 3e4); ax.legend(loc="upper center", fontsize=9, ncol=1, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Each cell switches suddenly, at its own time: all-or-none induction (Vilar et al. 2003, Fig. 3)",
                 x=0.01, ha="left", fontsize=13)
    fig.tight_layout()
    if save:
        fig.savefig(OUT / "fig_a4_vilar3.png", dpi=180)
    return fig


def fig_bimodality(plt=None, save=True):
    plt = plt or _style()
    R = json.loads((OUT / "vilar_results.json").read_text())
    trs = {"wt": Trajectory.from_parquet(OUT / "trajectory_vilar3a_sub.parquet"),
           "rho5": Trajectory.from_parquet(OUT / "trajectory_rho5_graded.parquet")}
    times = [0, 240, 480, 720, 1200, 1800]
    fig, axs = plt.subplots(2, len(times), figsize=(17, 5.6), sharex=True)
    bins = np.linspace(0, 3.4, 45)
    for r, (k, lab, col) in enumerate([("wt", f"Normal switch (ρ = 170), {R['T_sub_uM']:g} µM TMG", BLUE),
                                       ("rho5", f"Decoy strain (ρ = 5), {R['rho5_control']['TMG_uM']:.1f} µM TMG", ORANGE)]):
        for c, t in enumerate(times):
            ax = axs[r, c]
            v = np.log10(trs[k].at(t).Y_count.values + 1)
            ax.hist(v, bins=bins, color=col, alpha=0.85)
            b = bimodality_of_counts(trs[k].at(t).Y_count.values)
            ax.set_title(f"{t / 60:g} h" + (f"\ndip p = {b['dip_p']:.3f}" if b["dip_p"] is not None else ""), fontsize=10)
            ax.set_yticks([])
            if c == 0:
                ax.set_ylabel(lab.split(",")[0], fontsize=10)
            if r == 1:
                ax.set_xlabel("log10(LacY + 1)")
    fig.suptitle("Two separate groups of cells (top) vs one broad group (bottom): bimodal, not just noisy",
                 x=0.01, ha="left", fontsize=13)
    fig.text(0.01, -0.01, "Dip test: p < 0.05 means the histogram has more than one peak. Bottom row: same average LacY "
             "at 30 h, but no switch (ρ below the limit of 9).", fontsize=9, color=INK2)
    fig.tight_layout()
    if save:
        fig.savefig(OUT / "fig_a4_bimodality.png", dpi=180, bbox_inches="tight")
    return fig


def fig_hysteresis(plt=None, save=True):
    plt = plt or _style()
    R = json.loads((OUT / "vilar_results.json").read_text())
    h = pd.read_csv(OUT / "ssa_hysteresis.csv")
    p = dict(rn.DEFAULT_PARAMS); p.update(__import__("yaml").safe_load(YAML.read_text())["parameters"])
    off, on = rn.thresholds_tmg(0.0, p)
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.axvspan(off, on, color=GRID, alpha=0.8, lw=0)
    for hist, col, lab in [("uninduced", BLUE, "cells start OFF"), ("induced", ORANGE, "cells start ON")]:
        d = h[h.history == hist]
        ax.semilogx(d.TMG_uM, 100 * d.fraction_ON_20h, "o-", color=col, label=lab, ms=5)
    s = R["ssa_hysteresis"]
    ax.text(np.sqrt(off * on), 50, "memory\nzone\n(ODE)", ha="center", color=INK2)
    ax.set_xlabel("TMG outside the cell (µM)"); ax.set_ylabel("% of cells ON after 20 h")
    ax.set_title(f"1,000 simulated cells per point. ON cells keep their memory; OFF cells drift ON slowly\n50% switch points: ON at "
                 f"{s['switch_on_T50_uM']:.0f} µM, OFF at {s['switch_off_T50_uM']:.1f} µM (ODE 30 / 4.1)",
                 loc="left", fontsize=11.5)
    ax.legend(loc="center right")
    fig.tight_layout()
    if save:
        fig.savefig(OUT / "fig_a4_ssa_hysteresis.png", dpi=180)
    return fig


def fig_molecular(plt=None, save=True):
    """Molecular -> cellular level in one switching cell (exact SSA, from a4_validate_modes.py)."""
    plt = plt or _style()
    f = OUT / "exact_T40_from_OFF.parquet"
    if not f.exists():
        return None
    tr = Trajectory.from_parquet(f)
    ts = first_switch_times(tr)
    c = int(np.argsort(np.where(np.isfinite(ts), ts, 1e9))[min(3, np.isfinite(ts).sum() - 1)])
    d = tr.data[tr.data.cell_id == c]
    h = d.time.values / 60
    fig, axs = plt.subplots(4, 1, figsize=(8, 7), sharex=True)
    axs[0].step(h, d.P_free_count, color=PURPLE, where="post"); axs[0].set_ylabel("promoter\nfree?"); axs[0].set_yticks([0, 1])
    axs[1].plot(h, d.M_count, color=GREEN, lw=1); axs[1].set_ylabel("lac mRNA")
    axs[2].plot(h, d.Y_count, color=BLUE, lw=1.6); axs[2].set_ylabel("LacY")
    axs[3].plot(h, d.X_count, color=ORANGE, lw=1.2); axs[3].set_ylabel("TMG inside")
    for ax in axs:
        ax.axvline(ts[c] / 60, color=INK2, ls=":", lw=1)
    axs[3].set_xlabel("hours (40 µM TMG, exact simulation, sampled every 5 min)")
    axs[0].set_title("One cell switching ON: repressor falls off → mRNA bursts →\nmore LacY → more TMG inside → promoter stays free",
                     loc="left", fontsize=11.5)
    fig.tight_layout()
    if save:
        fig.savefig(OUT / "fig_a4_molecular.png", dpi=180)
    return fig


def postprocess():
    """Recompute every bimodality / summary number from the saved trajectories (de-quantised dip test)."""
    R = json.loads((OUT / "vilar_results.json").read_text())
    times = [0, 240, 480, 720, 1200, 1800]
    tr = Trajectory.from_parquet(OUT / "trajectory_vilar3a_sub.parquet")
    R["bimodality_wt"] = {str(t): bimodality_of_counts(tr.at(t).Y_count.values) |
                          {"fraction_ON": float((tr.at(t).Y_count.values > tr.meta["on_threshold_Y"]).mean())} for t in times}
    tr5 = Trajectory.from_parquet(OUT / "trajectory_rho5_graded.parquet")
    R["rho5_control"]["bimodality_rho5"] = {str(t): bimodality_of_counts(tr5.at(t).Y_count.values) for t in times}
    for key in R["vilar"]:
        t = Trajectory.from_parquet(OUT / f"trajectory_{key}.parquet")
        R["vilar"][key]["summary_stats"] = t.write_summary(OUT / f"summary_stats_{key}.json")
        if key == "vilar3a_sub":
            (ROOT / "outputs" / "summary_stats.json").write_text(json.dumps(R["vilar"][key]["summary_stats"], indent=2, default=float))
    t5 = tr5.write_summary(OUT / "summary_stats_rho5_graded.json")
    R["rho5_control"]["summary_stats"] = t5
    (OUT / "vilar_results.json").write_text(json.dumps(R, indent=2, default=float))
    return R


def figures():
    plt = _style()
    for f in (fig_vilar3, fig_bimodality, fig_hysteresis, fig_molecular):
        fig = f(plt)
        if fig is not None:
            plt.close(fig)


if __name__ == "__main__":
    if "--figures-only" not in sys.argv:
        R = simulate()
        R = postprocess()
        print(json.dumps({k: R[k] for k in ("ssa_hysteresis", "T_sub_uM")}, indent=1, default=float))
    figures()
