"""
A3 driver: fit -> hysteresis check -> figures -> circuit_model.yaml / SBML -> report numbers.
Usage: python src/run_a3.py <precise1k data root> <outputs dir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).parent))
import fit_parameters as fp  # noqa: E402
import reaction_network as rn  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 2,
    "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
})


def sweep(p, G, Ts, rho=None):
    q = dict(p, rho=rho) if rho else p
    up = np.array([rn.simulate(T, G, q, induced=False) for T in Ts])
    down = np.array([rn.simulate(T, G, q, induced=True) for T in Ts])
    return up, down


def switch_points(Ts, up, down, p, G):
    """Ozbudak's definition applied to the ODE: TMG where a history's state has jumped."""
    mid = np.sqrt(rn.y_max(p) * rn.c_crp(G, p) * rn.y_max(p) * rn.c_crp(G, p) / rn.rho_eff(G, p))
    on = Ts[np.argmax(up > mid)] if (up > mid).any() else np.nan
    off = Ts[np.where(down < mid)[0].max()] if (down < mid).any() else np.nan
    return off, on


def hill_coefficient(Ts, y):
    lo, hi = y.min(), y.max()
    f = (y - lo) / (hi - lo)
    t10, t90 = np.interp(0.1, f, Ts), np.interp(0.9, f, Ts)
    return np.log(81) / np.log(t90 / t10)


def lac_floor_analysis(root: Path, gat_ref, out: Path):
    """Why A2 cannot validate CRP -> lac: cAMP titration (K-12 'raman') vs IPTG-induced cells."""
    import pandas as pd
    from scipy import stats
    md = pd.read_csv(root / "k12_modulome" / "metadata_qc.csv", index_col=0)
    keep = md.project.isin(["raman", "srna_oe"])
    cnt = pd.read_csv(root / "k12_modulome" / "counts.csv", index_col=0)
    cnt = cnt[[x for x in md[keep].index if x in cnt.columns]]
    gi = pd.read_csv(root / "annotation" / "gene_info.csv", index_col=0)
    L = (gi.end - gi.start).abs().reindex(cnt.index) / 1000
    ok = L.notna() & (L > 0)
    rpk = cnt[ok].div(L[ok], axis=0); tpm = rpk / rpk.sum() * 1e6
    df = pd.DataFrame({"lacY": tpm.loc["b0343"], "gat": tpm.loc[["b2096", "b2095", "b2094"]].mean()})
    df = df.join(md[["project", "condition", "Supplement"]])
    ind = df[(df.project == "srna_oe") & df.Supplement.astype(str).str.contains("IPTG")]
    ram = df[df.project == "raman"]
    ratio_ind = float((ind.lacY / ind.gat).median())
    lr = stats.linregress(np.log2(ind.gat), np.log2(ind.lacY))
    order = ["del_cyaA", "del_cyaA_01mM", "del_cyaA_05mM", "del_cyaA_1mM", "wt"]
    labels = ["cyaA KO\n0 cAMP", "KO +0.1\nmM cAMP", "KO +0.5\nmM", "KO +1\nmM", "Wild\ntype"]
    g = ram.groupby("condition")[["gat", "lacY"]].median().loc[order]
    pred = ratio_ind / 170 * g.gat
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    x = np.arange(len(order))
    ax.semilogy(x, g.gat / gat_ref, "o-", color=AQUA, label="CRP activity (gat, ÷ reference)")
    ax.semilogy(x, g.lacY, "o-", color=BLUE, label="lacY measured (TPM)")
    ax.semilogy(x, pred, "o--", color=BLUE, alpha=0.45, label="lacY expected if OFF lac ∝ CRP (TPM)")
    ax.axhspan(0.005, 1.0, color=GRID, alpha=0.8, lw=0)
    ax.text(1.55, 0.6, "below RNA-seq detection floor (~1 TPM)", color=INK2, fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0.005, 10)
    ax.set_title("CRP rises 40×, but OFF-state lac stays under the floor", loc="left", fontsize=11.5)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout(); fig.savefig(out / "fig_a3_lac_floor.png", dpi=200); plt.close(fig)
    return {"induced_lacY_vs_gat_loglog_slope": lr.slope, "induced_slope_se": lr.stderr, "induced_r": lr.rvalue,
            "induced_n": int(len(ind)),
            "raman_median_by_condition": g.round(2).to_dict(), "raman_expected_basal_lacY": pred.round(3).to_dict()}


def fig_outputs(otab, oval, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    lim = [-6.5, 2]
    ax.plot(lim, lim, color=GRID, lw=1)
    weak = otab.R2_all_samples < 0.06
    ax.scatter(otab.a2_ko_effect_log2[~weak], otab.pred_ko_effect_log2[~weak], s=40, color=BLUE, zorder=3)
    ax.scatter(otab.a2_ko_effect_log2[weak], otab.pred_ko_effect_log2[weak], s=40, facecolor="white",
               edgecolor=BLUE, lw=1.5, zorder=3, label="weak fit (R² < 0.06)")
    for r in otab.itertuples():
        lab = {"malK": "malK, lamB", "flhD": "flhC, flhD", "rbsA": "rbsA", "rbsB": "rbsB", "rbsC": "rbsC",
               "uhpT": "uhpT", "gyrA": "gyrA", "fecE": "fecE"}.get(r.target)
        if lab:
            ax.annotate(lab, (r.a2_ko_effect_log2, r.pred_ko_effect_log2), xytext=(4, 3),
                        textcoords="offset points", fontsize=8, color=INK2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("A2 measured: crp knockout effect (log2)")
    ax.set_ylabel("Model predicted (log2)")
    ax.set_title("21 CRP outputs predict A2's knockout effects", loc="left", fontsize=11.5)
    ax.legend(loc="upper left", fontsize=9)
    ax = axes[1]
    lim = [-7, 2]
    ax.plot(lim, lim, color=GRID, lw=1)
    cmap = {"del_cyaA": INK, "del_cyaA_01mM": ORANGE, "del_cyaA_05mM": BLUE, "del_cyaA_1mM": AQUA}
    lab = {"del_cyaA": "cyaA KO, 0 cAMP", "del_cyaA_01mM": "+0.1 mM", "del_cyaA_05mM": "+0.5 mM", "del_cyaA_1mM": "+1 mM"}
    for cnd, col in cmap.items():
        v = oval[oval.condition == cnd]
        ax.scatter(v.obs_log2FC_vs_wt, v.pred_log2FC_vs_wt, s=26, color=col, label=lab[cnd], zorder=3,
                   edgecolor="white", lw=0.8)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Measured change vs wild type (log2)")
    ax.set_ylabel("Model predicted (log2)")
    ax.set_title("Unseen test: cAMP titration in a cyaA knockout", loc="left", fontsize=11.5)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout(); fig.savefig(out / "fig_a3_crp_outputs.png", dpi=200); plt.close(fig)


def main(data_root, out):
    out = Path(out)
    p, fitted, tab, rho_tab, p1k, k12 = fp.run_all(Path(data_root), out)
    res = {"rho_by_condition": fitted["derived"]["rho_by_condition"]}

    # ---------- 0. A2-validated CRP outputs + why lac cannot be A2-validated
    import crp_outputs as co
    otab, oval, osum, curves = co.run(Path(data_root), out)
    outputs = {r.target: {"basal_TPM": r.basal_TPM, "amplitude_TPM": r.amplitude_TPM, "K_crp": r.K_crp,
                          "hill_n": r.hill_n, "direction": r.direction,
                          "R2_all_samples": r.R2_all_samples, "a2_ko_effect_log2": r.a2_ko_effect_log2,
                          "pred_ko_effect_log2": r.pred_ko_effect_log2,
                          "fit_quality": "weak" if r.R2_all_samples < 0.06 else "ok"}
               for r in otab.itertuples()}
    res["crp_outputs"] = osum
    lac_floor = lac_floor_analysis(Path(data_root), osum["gat_reference_TPM"], out)
    res["lac_floor"] = lac_floor
    fig_outputs(otab, oval, out)

    # ---------- 1. dynamic hysteresis check (20 h ODE from each history)
    Ts = np.logspace(-0.3, 3.3, 90)
    sweeps = {}
    for G in [0.0, 1000.0]:
        up, down = sweep(p, G, Ts)
        off, on = switch_points(Ts, up, down, p, G)
        an_off, an_on = rn.thresholds_tmg(G, p)
        sweeps[G] = (up, down)
        res[f"ode_sweep_G{int(G)}"] = {"switch_off_ode": off, "switch_on_ode": on,
                                       "switch_off_analytic": an_off, "switch_on_analytic": an_on,
                                       "on_off_fold_expression": float(down.max() / up.min())}
    i = np.argmin(abs(Ts - 15))
    res["fold_induced_vs_uninduced_at_15uM_G0"] = float(sweeps[0.0][1][i] / sweeps[0.0][0][i])

    # ---------- 2. cusp: smallest rho with bistability
    rhos = np.linspace(5, 15, 1001)
    bist = [rn.switching_drives(r, p["n"]) is not None for r in rhos]
    res["cusp_rho"] = float(rhos[np.argmax(bist)])

    # ---------- 3. titration strains at 1 mM glucose (Ozbudak Fig 4a-c)
    strains = {}
    Tg = np.logspace(0, 4, 120)
    for rho in [170, 50, 5]:
        up, down = sweep(p, 1000.0, Tg, rho=rho)
        hyst = np.max(np.abs(np.log(down / up))) > 0.5
        strains[rho] = (up, down)
        res[f"strain_rho{rho}"] = {"hysteretic": bool(hyst),
                                    "hill_coefficient_up": float(hill_coefficient(Tg, up))}

    # ---------- figures
    figdir = out
    # (a) hysteresis at 0 glucose
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    up, down = sweeps[0.0]
    off, on = rn.thresholds_tmg(0.0, p)
    ax.axvspan(off, on, color=GRID, alpha=0.7, lw=0)
    ax.loglog(Ts, up, color=BLUE, label="Cells start OFF")
    ax.loglog(Ts, down, color=ORANGE, ls=(0, (4, 2)), label="Cells start ON")
    for x in (3, 30):
        ax.axvline(x, color=INK2, lw=0.8, ls=":")
        ax.text(x * 1.06, 40, f"{x} µM\n(paper)", color=INK2, fontsize=9, va="center")
    ax.text(np.sqrt(off * on), 350, "memory\nzone", ha="center", color=INK2, fontsize=10)
    ax.set_xlabel("TMG outside the cell (µM)")
    ax.set_ylabel("LacY per cell (molecules)")
    ax.set_title("A3 model shows memory: two answers between 4 and 30 µM", loc="left", fontsize=11.5)
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(figdir / "fig_a3_hysteresis.png", dpi=200); plt.close(fig)

    # (b) phase diagram vs Ozbudak Fig 2c
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    Gs = np.logspace(-1, 3.3, 200)
    th = np.array([rn.thresholds_tmg(G, p) for G in Gs])
    Gplot = np.r_[0.03, Gs]
    th = np.vstack([rn.thresholds_tmg(0.0, p), th])
    ax.fill_between(Gplot, th[:, 0], th[:, 1], color=GRID, alpha=0.9, lw=0, label="Model: bistable")
    ax.plot(Gplot, th[:, 1], color=BLUE, lw=1.5)
    ax.plot(Gplot, th[:, 0], color=ORANGE, lw=1.5)
    d = fp.FIG2C.astype(float); d.loc[d.glucose_uM == 0, "glucose_uM"] = 0.03
    ax.scatter(d.glucose_uM, d.switch_on_uM, s=46, color=BLUE, zorder=3, edgecolor="white", lw=1.5,
               label="Paper switch-ON (used to fit)")
    ax.scatter(d.glucose_uM, d.switch_off_uM, s=46, facecolor="white", edgecolor=ORANGE, lw=1.8, zorder=3,
               label="Paper switch-OFF (held out)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([0.03, 0.1, 1, 10, 100, 1000]); ax.set_xticklabels(["0", "0.1", "1", "10", "100", "1000"])
    ax.set_xlabel("Glucose outside the cell (µM)")
    ax.set_ylabel("TMG outside the cell (µM)")
    ax.set_title("Model reproduces Ozbudak's phase diagram (Fig. 2c)", loc="left", fontsize=11.5)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout(); fig.savefig(figdir / "fig_a3_phase_diagram.png", dpi=200); plt.close(fig)

    # (c) rho: where every estimate sits relative to the cusp
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    rows = [(f"{r.gene} RNA, vs {'same-batch WT' if 'same' in r.uninduced_reference else 'PRECISE-1K WT'}",
             r.rho, r.ci_lo, r.ci_hi) for r in rho_tab.itertuples()]
    rows += [("Fig. 2c thresholds (free fit)", fitted["rho_free_fit_to_fig2c"], np.nan, np.nan),
             ("Ozbudak 2004 (GFP protein)", 170, 140, 200)]
    if p.get("rho_crp_coupling"):
        rows += [("Model, crp knockout (Kuhlman looping)", rn.rho_eff(0.0, p, crp_knockout=True), np.nan, np.nan)]
    for k, (lab, v, lo, hi) in enumerate(rows):
        c = INK if "Ozbudak" in lab or "Fig" in lab else (AQUA if "Model" in lab else BLUE)
        if not np.isnan(lo):
            ax.plot([lo, hi], [k, k], color=c, lw=1.5)
        ax.scatter(v, k, color=c, s=40, zorder=3)
    ax.axvspan(1, 9, color=ORANGE, alpha=0.12, lw=0)
    ax.axvline(9, color=ORANGE, lw=1)
    ax.text(8.5, -0.3, "no switch\n(ρ < 9)", ha="right", va="bottom", color=ORANGE, fontsize=9)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xscale("log"); ax.set_xlim(3, 1000)
    ax.set_xlabel("Repression factor ρ (fully ON ÷ OFF)")
    ax.set_title("Every estimate of ρ is far above the switch limit of 9", loc="left", fontsize=11.5)
    ax.grid(axis="y", visible=False)
    fig.tight_layout(); fig.savefig(figdir / "fig_a3_rho.png", dpi=200); plt.close(fig)

    # (d) titration strains
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, rho in zip(axes, [170, 50, 5]):
        up, down = strains[rho]
        ax.loglog(Tg, up, color=BLUE, label="start OFF")
        ax.loglog(Tg, down, color=ORANGE, ls=(0, (4, 2)), label="start ON")
        kind = "switch (memory)" if res[f"strain_rho{rho}"]["hysteretic"] else "dimmer (no memory)"
        ax.set_title(f"ρ = {rho}: {kind}", loc="left", fontsize=11)
        ax.set_xlabel("TMG (µM), 1 mM glucose")
    axes[0].set_ylabel("LacY per cell"); axes[0].legend(loc="upper left", fontsize=9)
    fig.tight_layout(); fig.savefig(figdir / "fig_a3_titration.png", dpi=200); plt.close(fig)

    # combined sanity-check figure required by the brief
    import matplotlib.image as mpimg
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))
    for ax, f in zip(axs.flat, ["fig_a3_hysteresis", "fig_a3_phase_diagram", "fig_a3_rho", "fig_a3_titration"]):
        ax.imshow(mpimg.imread(figdir / f"{f}.png")); ax.axis("off")
    fig.tight_layout(); fig.savefig(figdir / "hysteresis_check.png", dpi=150); plt.close(fig)

    # ---------- circuit_model.yaml
    q = {k: (v if isinstance(v, str) or v is None else float(v)) for k, v in p.items()}
    model = {
        "name": "crp_lac_circuit_A3",
        "units": {"time": "min", "species": "molecules per cell", "TMG_ext": "uM", "glucose": "uM"},
        "inputs": {"TMG_ext": "extracellular TMG (uM)", "glucose": "extracellular glucose (uM)"},
        "edges": [
            {"source": "CRP-cAMP", "target": "P_lac", "sign": "+",
             "provenance": "published_intervention: Kuhlman et al. 2007 PNAS 104:6043 Fig 1b (dcyaA dcpdA dlacY, "
                           "cAMP titration at 1 mM IPTG; 201x rise 0 -> 1 mM, flat in dcrp); also RegulonDB + Ozbudak Fig 3b",
             "evidence": "published_intervention", "a2_validated": False,
             "a2_note": "lacZYA showed no CRP-perturbation effect in PRECISE-1K (no inducer in any sample)",
             "quantitative_check": "output close to linear in CRP occupancy (Hill m = 1.15, CI 1.03-1.28; "
                                   "paper 1.0 +/- 0.1), as in R3"},
            {"source": "CRP-cAMP", "target": "LacI-P_lac loop", "sign": "+ (tighter repression)",
             "provenance": "published_intervention: Kuhlman et al. 2007 Fig 2a (looping L vs cAMP; Omega = 12.4, "
                           "CI 10.2-15.4 refit, paper 10.3)",
             "evidence": "published_intervention", "a2_validated": False,
             "model_term": "rho_eff in R1; default applies it to crp knockout only (rho 170 -> ~22); "
                           "glucose form (rho_crp_coupling: glucose) fails Ozbudak Fig 2c, see fitted_params sensitivity"},
            {"source": "LacI", "target": "P_lac", "sign": "-", "provenance": "prior: RegulonDB + Ozbudak Fig 3c",
             "a2_validated": False},
            {"source": "TMG_in", "target": "LacI", "sign": "-", "provenance": "prior: literature", "a2_validated": False},
            {"source": "LacY", "target": "TMG_in", "sign": "+", "provenance": "prior: Ozbudak Fig 3d", "a2_validated": False},
            {"source": "glucose", "target": "CRP-cAMP", "sign": "-", "provenance": "PRECISE-1K (gat 2.3x lower on glucose) + Ozbudak Fig 3a",
             "a2_validated": False},
            {"source": "glucose", "target": "LacY", "sign": "-", "provenance": "prior: inducer exclusion, Ozbudak Fig 3e",
             "a2_validated": False},
        ],
        "species": {"P_free": {"initial": 0}, "P_rep": {"initial": 1}, "M": {"initial": 0},
                    "Y": {"initial": round(rn.y_max(p) / rn.rho_eff(0.0, p))}, "X": {"initial": 0}},
        "parameters": q,
        "derived": fitted["derived"],
        "input_functions": {
            "c_crp": "crp_min + (1 - crp_min) / (1 + glucose / crp_K)",
            "beta_G": "excl_min + (1 - excl_min) / (1 + glucose / excl_K)",
            "beta_T": "TMG_ext ** uptake_h   (power law, uM units)",
            "rho_eff": "1 + (rho - 1) * (1 + L(theta)) / (1 + L(crp_occ_ref)), L(theta) = loop_L0 * (1 + (loop_omega - 1) * theta); "
                       "theta = 0 for crp knockout; with rho_crp_coupling 'glucose' also theta = crp_occ_ref * c_crp(glucose); "
                       "otherwise rho_eff = rho",
        },
        "reactions": [
            {"id": "R1_repress", "reactants": {"P_free": 1}, "products": {"P_rep": 1},
             "propensity": "k_rep_max(glucose) / (1 + (X/x0)**n) * P_free,  k_rep_max = (rho_eff(glucose) - 1) * k_unrep",
             "k_rep_max_no_glucose": rn.k_rep_max(p, 0.0), "k_rep_max_glucose_1mM": rn.k_rep_max(p, 1000.0)},
            {"id": "R2_derepress", "reactants": {"P_rep": 1}, "products": {"P_free": 1}, "propensity": "k_unrep * P_rep"},
            {"id": "R3_transcribe", "reactants": {"P_free": 1}, "products": {"P_free": 1, "M": 1},
             "propensity": "k_tx * c_crp(glucose) * P_free"},
            {"id": "R4_mrna_decay", "reactants": {"M": 1}, "products": {}, "propensity": "d_m * M"},
            {"id": "R5_translate", "reactants": {"M": 1}, "products": {"M": 1, "Y": 1}, "propensity": "k_tl * M"},
            {"id": "R6_lacY_dilute", "reactants": {"Y": 1}, "products": {}, "propensity": "ln(2)/doubling_min * Y"},
            {"id": "R7_uptake", "reactants": {"Y": 1}, "products": {"Y": 1, "X": 1},
             "propensity": "k_up * beta_T(TMG_ext) * beta_G(glucose) * Y"},
            {"id": "R8_tmg_loss", "reactants": {"X": 1}, "products": {}, "propensity": "d_x * X"},
        ],
        "validation": {"ozbudak_fig2c_switch_on_rmse_log10": fitted["fit"]["switch_on_rmse_log10"],
                       "ozbudak_fig2c_switch_off_rmse_log10_heldout": fitted["fit"]["switch_off_rmse_log10_heldout"],
                       "cusp_rho": res["cusp_rho"],
                       "rho_by_condition": res["rho_by_condition"],
                       "crp_outputs_ko_effect_r": osum["ko_effect_pearson_r"],
                       "crp_outputs_raman_r": osum["raman_pearson_r"]},
        "crp_outputs": {
            "description": "A2-validated CRP targets as outputs; mRNA (TPM) = basal + amplitude * Hill(c_crp(glucose)); "
                           "repressed targets use 1 - Hill. Same c_crp drives the lac promoter; crp knockout -> c = 0. "
                           "Fitted on PRECISE-1K (n = 1009, DGF-298 excluded); c read from gat TPM / reference.",
            "gat_reference_TPM": osum["gat_reference_TPM"],
            "tpm_to_molecules": 0.0024,
            "targets": outputs,
        },
    }
    model["edges"] += [{"source": "CRP-cAMP", "target": k, "sign": "+" if v["direction"] == "activated" else "-",
                        "provenance": "A2 causal-mediation test (PRECISE-1K crp project)", "a2_validated": True,
                        "fit_quality": v["fit_quality"]} for k, v in outputs.items()]
    model = json.loads(json.dumps(model, default=float))
    (out / "circuit_model.yaml").write_text(yaml.safe_dump(model, sort_keys=False))

    # ---------- SBML (COPASI / GillesPy2) at a reference condition inside the bistable band
    T_ref, G_ref = 15.0, 0.0
    gm = rn.to_gillespy2(p, T_ref, G_ref, outputs=outputs)
    try:
        import gillespy2
        gillespy2.export_SBML(gm, filename=str(out / "circuit_model_TMG15_G0.xml"))
        res["sbml_export"] = "ok"
    except Exception as e:  # noqa: BLE001
        res["sbml_export"] = f"failed: {e}"

    # ---------- GillesPy2 smoke test (A4 preview, not the A4 result)
    try:
        gm.timespan(np.linspace(0, 600, 61))
        try:
            r, solver = gm.run(number_of_trajectories=1, seed=7), "default (C++ SSA)"
        except Exception:  # noqa: BLE001  C++ solver needs a working build toolchain; fall back to pure Python
            r, solver = gm.run(solver=gillespy2.NumPySSASolver, number_of_trajectories=1, seed=7), "NumPySSASolver"
        res["gillespy2_smoke_test"] = {"ok": True, "final_LacY": int(r[0]["Y"][-1]), "solver": solver}
    except Exception as e:  # noqa: BLE001
        res["gillespy2_smoke_test"] = {"ok": False, "error": str(e)[:300]}

    (out / "a3_checks.json").write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
