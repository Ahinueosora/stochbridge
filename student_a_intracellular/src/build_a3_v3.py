"""
A3 v3 build: the final parameterized CRP/lac circuit on the v3 network (A2 + RegulonDB 14.5.1 + Kuhlman 2007 + Ozbudak 2004).

Reuses the A3 fit (outputs/a3/fitted_params.json, outputs/a3/crp_outputs_fit.csv; refit needs PRECISE-1K: run_a3.py) and writes
    outputs/a3/circuit_model.yaml        full spec: species, reactions, rates, initial conditions, edges + evidence
    outputs/a3/circuit_model.xml         SBML L3V2 of the same model (COPASI + GillesPy2), TMG 15 uM, 0 glucose, cells start OFF
    outputs/a3/validated_edges.json      every circuit edge with where its validation comes from
    outputs/a3/hysteresis_check.png      (+ the four panel PNGs)
    outputs/a3/a3_checks_v3.json         all numbers quoted in parameter_fit_report.md

Usage: python src/build_a3_v3.py <network_v3_edges.csv> [outputs dir, default outputs/a3]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import reaction_network as rn  # noqa: E402
from fit_parameters import FIG2C  # noqa: E402
from run_a3 import sweep, switch_points, hill_coefficient, BLUE, ORANGE, AQUA, INK, INK2, GRID  # noqa: E402

RED = "#c4314b"
A2_GENES = ["hflC", "tdcB", "uhpT", "cas3", "rbsC", "yahN", "fecE", "mhpF", "flhC", "flhD", "araJ", "rbsA", "nadC",
            "galK", "rbsB", "gyrA", "rbsR", "lamB", "malK", "preA", "idnO"]
SPECIES_DOC = {
    "P_free": "lac promoter not bound by LacI (0 or 1 copy)",
    "P_rep": "lac promoter bound by LacI (0 or 1 copy); P_free + P_rep = 1",
    "M": "lacZYA mRNA (molecules)",
    "Y": "LacY permease (molecules)",
    "X": "internal TMG (molecules); u = X/x0 inactivates LacI",
}


# ------------------------------------------------------------------ helpers
def load_inputs(out: Path):
    fitted = json.loads((out / "fitted_params.json").read_text())
    p = dict(rn.DEFAULT_PARAMS)
    p.update(fitted["params"])
    p.setdefault("crp_free_basal", rn.DEFAULT_PARAMS["crp_free_basal"])
    ofit = pd.read_csv(out / "crp_outputs_fit.csv")
    outputs = {r.target: {"basal_TPM": float(r.basal_TPM), "amplitude_TPM": float(r.amplitude_TPM),
                          "K_crp": float(r.K_crp), "hill_n": float(r.hill_n), "direction": r.direction,
                          "R2_all_samples": float(r.R2_all_samples), "a2_ko_effect_log2": float(r.a2_ko_effect_log2),
                          "pred_ko_effect_log2": float(r.pred_ko_effect_log2),
                          "fit_quality": "weak" if r.R2_all_samples < 0.06 else "ok"}
               for r in ofit.itertuples()}
    return p, fitted, outputs


def rate_constants(p, T, G, crp_knockout=False):
    """Numbers that go into the propensities at one condition (inputs fixed for a run)."""
    return {"k_rep_max": rn.k_rep_max(p, G, crp_knockout), "k_unrep": p["k_unrep"], "x0": p["x0"], "n": p["n"],
            "k_tx_eff": p["k_tx"] * rn.c_eff(G, p, crp_knockout), "d_m": p["d_m"], "k_tl": p["k_tl"],
            "d_y": rn.d_y(p), "k_up_eff": p["k_up"] * float(rn.beta_T(T, p)) * rn.beta_G(G, p), "d_x": p["d_x"]}


def initial_counts(p, G, induced, crp_knockout=False):
    P_free, M, Y, X = rn.initial_state(p, induced, G, crp_knockout)
    if induced:
        return {"P_free": 1, "P_rep": 0, "M": int(round(M)), "Y": int(round(Y)), "X": int(round(X))}
    return {"P_free": 0, "P_rep": 1, "M": 0, "Y": int(round(Y)), "X": 0}


def output_rates(p, G, outputs, tpm_to_molecules=0.0024, crp_knockout=False):
    tpm = rn.crp_output_tpm(G, outputs, p, crp_knockout=crp_knockout)
    return {g: float(v * tpm_to_molecules * p["d_m"]) for g, v in tpm.items()}


# ------------------------------------------------------------------ SBML (COPASI + GillesPy2)
REACTIONS = [  # id, reactants, products, kinetic law (molecules / min)
    ("R1_repress", {"P_free": 1}, {"P_rep": 1}, "k_rep_max / (1 + (X / x0)^n) * P_free"),
    ("R2_derepress", {"P_rep": 1}, {"P_free": 1}, "k_unrep * P_rep"),
    ("R3_transcribe", {"P_free": 1}, {"P_free": 1, "M": 1}, "k_tx_eff * P_free"),
    ("R4_mrna_decay", {"M": 1}, {}, "d_m * M"),
    ("R5_translate", {"M": 1}, {"M": 1, "Y": 1}, "k_tl * M"),
    ("R6_lacY_dilute", {"Y": 1}, {}, "d_y * Y"),
    ("R7_uptake", {"Y": 1}, {"Y": 1, "X": 1}, "k_up_eff * Y"),
    ("R8_tmg_loss", {"X": 1}, {}, "d_x * X"),
]


def write_sbml(path: Path, p, T, G, outputs=None, induced=False, crp_knockout=False):
    import libsbml
    doc = libsbml.SBMLDocument(3, 2)
    m = doc.createModel()
    m.setId("crp_lac_circuit_v3")
    m.setName(f"CRP/lac circuit A3 v3, TMG {T:g} uM, glucose {G:g} uM, "
              f"{'induced' if induced else 'uninduced'} start{', crp knockout' if crp_knockout else ''}")
    # units: molecules (item) per cell, minutes
    def unit(uid, kinds):
        ud = m.createUnitDefinition(); ud.setId(uid)
        for kind, exp, scale, mult in kinds:
            u = ud.createUnit(); u.setKind(kind); u.setExponent(exp); u.setScale(scale); u.setMultiplier(mult)
    unit("minute", [(libsbml.UNIT_KIND_SECOND, 1, 0, 60)])
    unit("per_minute", [(libsbml.UNIT_KIND_SECOND, -1, 0, 60)])
    unit("item_per_minute", [(libsbml.UNIT_KIND_ITEM, 1, 0, 1), (libsbml.UNIT_KIND_SECOND, -1, 0, 60)])
    m.setTimeUnits("minute"); m.setSubstanceUnits("item"); m.setExtentUnits("item"); m.setVolumeUnits("litre")
    c = m.createCompartment(); c.setId("cell"); c.setSize(1.0); c.setConstant(True); c.setSpatialDimensions(3)

    init = initial_counts(p, G, induced, crp_knockout)
    species = dict(init)
    orates = output_rates(p, G, outputs, crp_knockout=crp_knockout) if outputs else {}
    for g in orates:
        species[f"mRNA_{g}"] = 0
    for sid, val in species.items():
        s = m.createSpecies(); s.setId(sid); s.setCompartment("cell"); s.setInitialAmount(val)
        s.setHasOnlySubstanceUnits(True); s.setBoundaryCondition(False); s.setConstant(False); s.setSubstanceUnits("item")

    pars = rate_constants(p, T, G, crp_knockout)
    pars.update({f"k_{g}": v for g, v in orates.items()})
    pars.update({"TMG_ext_uM": T, "glucose_uM": G})          # documentation: the condition these rates were computed at
    for pid, val in pars.items():
        q = m.createParameter(); q.setId(pid); q.setValue(float(val)); q.setConstant(True)
    rxns = list(REACTIONS)
    for g in orates:
        rxns += [(f"tx_{g}", {}, {f"mRNA_{g}": 1}, f"k_{g}"), (f"deg_{g}", {f"mRNA_{g}": 1}, {}, f"d_m * mRNA_{g}")]
    for rid, reac, prod, law in rxns:
        r = m.createReaction(); r.setId(rid); r.setReversible(False)
        for sid, st in reac.items():
            sr = r.createReactant(); sr.setSpecies(sid); sr.setStoichiometry(st); sr.setConstant(True)
        for sid, st in prod.items():
            sr = r.createProduct(); sr.setSpecies(sid); sr.setStoichiometry(st); sr.setConstant(True)
        import re
        for sid in sorted(set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", law)) & set(species) - set(reac) - set(prod)):
            mr = r.createModifier(); mr.setSpecies(sid)          # e.g. X in R1 (inducer inactivates LacI)
        kl = r.createKineticLaw(); kl.setMath(libsbml.parseL3Formula(law))
    libsbml.writeSBMLToFile(doc, str(path))
    return doc


def check_sbml(path: Path, p, T, G):
    """Validate with libSBML, load + run in GillesPy2 (SSA) and COPASI (deterministic + stochastic);
    compare COPASI's deterministic run with the Python ODE."""
    import libsbml
    res = {}
    doc = libsbml.readSBMLFromFile(str(path))
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
    n_err = doc.checkConsistency()
    errs = [doc.getError(i) for i in range(n_err)]
    res["libsbml_errors"] = sum(e.getSeverity() >= libsbml.LIBSBML_SEV_ERROR for e in errs)
    res["libsbml_warnings"] = sum(e.getSeverity() == libsbml.LIBSBML_SEV_WARNING for e in errs)

    import gillespy2
    gm, gerr = gillespy2.import_SBML(str(path))
    for sp in gm.listOfSpecies.values():        # import_SBML makes species continuous; SSA needs discrete counts
        sp.mode = "discrete"
        sp.initial_value = int(round(sp.initial_value))
    gm.timespan(gillespy2.TimeSpan.linspace(t=1200, num_points=121))
    r = gm.run(solver=gillespy2.SSACSolver, number_of_trajectories=20, seed=11)
    fin = np.array([t["Y"][-1] for t in r])
    res["gillespy2"] = {"import_messages": len(gerr), "species": len(gm.listOfSpecies),
                        "reactions": len(gm.listOfReactions), "solver": "SSACSolver (discrete)",
                        "ssa_20_cells_LacY_at_20h": {"min": int(fin.min()), "median": float(np.median(fin)),
                                                     "max": int(fin.max()), "n_ON_over_300": int((fin > 300).sum())}}

    import basico
    basico.load_model(str(path))
    det = basico.run_time_course(duration=1200, intervals=120, method="deterministic", use_numbers=True)
    sto_fin = [basico.run_time_course(duration=1200, intervals=12, method="directMethod", use_numbers=True,
                                      seed=s, use_seed=True)["Y"].values[-1] for s in range(1, 11)]
    from scipy.integrate import solve_ivp
    y0 = rn.initial_state(p, False, G)
    py = solve_ivp(rn.rhs, (0, 1200), y0, args=(T, G, p), method="LSODA", rtol=1e-8, atol=1e-10,
                   t_eval=det.index.values)
    copasi_y = det["Y"].values
    res["copasi"] = {"loaded": True, "species": len(basico.get_species()), "reactions": len(basico.get_reactions()),
                     "deterministic_LacY_at_20h": float(copasi_y[-1]), "python_ode_LacY_at_20h": float(py.y[2, -1]),
                     "max_rel_diff_LacY": float(np.max(np.abs(copasi_y - py.y[2]) / np.maximum(py.y[2], 1))),
                     "stochastic_10_cells_LacY_at_20h": [int(v) for v in sto_fin],
                     "units": basico.get_model_units()["quantity_unit"]}
    basico.remove_datamodel(basico.get_current_model())
    return res


# ------------------------------------------------------------------ edges + validation sources
def edge_tables(net_csv: Path, outputs, ext_json: Path | None = None):
    E = pd.read_csv(net_csv).fillna("")
    ext = {}
    if ext_json and Path(ext_json).exists():       # RegulonDB grades per CRP target (validated_edges_external_v2.json)
        ext = {(e["regulator"], e["target"]): e for e in json.loads(Path(ext_json).read_text())}

    def tags(ev, tier):
        t = []
        for e in ev.split(";"):
            e = e.strip()
            if e == "A2_causal": t.append("A2")
            elif e == "Kuhlman2007_intervention": t.append("Kuhlman2007_intervention")
            elif e.startswith("RegulonDB_Confirmed_grade1"): t.append("RegulonDB_Confirmed_grade1")
            elif e.startswith("RegulonDB_"): t.append(e)
            elif e == "Ozbudak2004": t.append("Ozbudak2004_measured")
        return sorted(set(t))

    core_map = {("CRP", "lac promoter (lacZYA)"): ("CRP-cAMP", "P_lac", "R3 (k_tx * c_crp)"),
                ("LacI", "lac promoter (lacZYA)"): ("LacI", "P_lac", "R1/R2 (rho)"),
                ("CRP", "LacI DNA loop"): ("CRP-cAMP", "LacI-P_lac loop", "rho_eff (crp knockout: rho 170 -> 22)"),
                ("inducer inside cell", "LacI"): ("TMG_in", "LacI", "R1 (1 + (X/x0)^n)"),
                ("LacY", "inducer inside cell"): ("LacY", "TMG_in", "R7 (uptake)"),
                ("inducer (TMG/IPTG)", "inducer inside cell"): ("TMG_ext", "TMG_in", "R7 (beta_T)"),
                ("lac promoter (lacZYA)", "LacY"): ("P_lac", "LacY", "R3-R5"),
                ("glucose", "cAMP"): ("glucose", "CRP-cAMP", "c_crp(glucose) in R3"),
                ("cAMP", "CRP"): ("cAMP", "CRP-cAMP", "folded into c_crp (Hill m = 1)"),
                ("glucose", "LacY"): ("glucose", "LacY activity", "beta_G in R7")}
    circuit, context = [], []
    for e in E.itertuples():
        k = (e.source, e.target)
        rec = {"source": e.source, "target": e.target, "sign": e.sign, "validated_by": tags(e.evidence, e.tier),
               "evidence": e.evidence.split(";") if e.evidence else [], "tier": e.tier,
               "strength": e.strength or None, "params": json.loads(e.params) if e.params else None,
               "note": e.note or None}
        if k in core_map:
            s, t, term = core_map[k]
            rec.update({"model_source": s, "model_target": t, "model_term": term, "simulated": True,
                        "a2_validated": False})
            circuit.append(rec)
        elif e.source == "CRP" and e.target in outputs:
            o = outputs[e.target]
            x = ext.get(("CRP", e.target))
            if x:
                rec["regulondb"] = {"effect": x["trn_effect"], "confidence": x["regulondb_confidence"],
                                    "causal_grade": x["causal_grade"], "pmids": x["pmids"],
                                    "sign_agrees_with_A2": x["agrees_with_a2"]}
                if x["externally_validated"]:
                    rec["validated_by"] = sorted(set(rec["validated_by"]) | {f"RegulonDB_grade{x['causal_grade']}"})
            rec.update({"model_source": "CRP-cAMP", "model_target": f"mRNA_{e.target}", "simulated": True,
                        "model_term": "Hill output of c_crp (tx_/deg_ reactions)", "a2_validated": True,
                        "a2_direction": o["direction"], "a2_ko_effect_log2": o["a2_ko_effect_log2"],
                        "fit_quality": o["fit_quality"], "sign_conflict_with_RegulonDB": "SIGN CONFLICT" in rec["note"] if rec["note"] else False})
            circuit.append(rec)
        else:
            rec.update({"simulated": False})
            context.append(rec)
    return circuit, context


# ------------------------------------------------------------------ figures
def figures(p, out: Path, res, sweeps, strains, Ts, Tg, ko_sweep):
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "lines.linewidth": 2, "legend.frameon": False, "figure.facecolor": "white"})
    # (a) hysteresis, wild type + crp knockout
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    up, down = sweeps[0.0]
    off, on = rn.thresholds_tmg(0.0, p)
    ax.axvspan(off, on, color=GRID, alpha=0.7, lw=0)
    ax.loglog(Ts, up, color=BLUE, label="wild type, cells start OFF")
    ax.loglog(Ts, down, color=ORANGE, ls=(0, (4, 2)), label="wild type, cells start ON")
    ax.loglog(Ts, ko_sweep[0], color=RED, lw=1.6, label="crp deleted (stays OFF)")
    for x in (3, 30):
        ax.axvline(x, color=INK2, lw=0.8, ls=":")
        ax.text(x * 1.06, 0.05, f"{x} µM\n(paper)", color=INK2, fontsize=9, va="center")
    ax.text(np.sqrt(off * on), 300, "memory\nzone", ha="center", color=INK2, fontsize=10)
    ax.set_ylim(0.01, 4000)
    ax.set_xlabel("TMG outside the cell (µM)"); ax.set_ylabel("LacY per cell (molecules)")
    ax.set_title("Memory between 4 and 30 µM; no switch without CRP", loc="left", fontsize=11.5)
    ax.legend(loc="center right", fontsize=8.5)
    fig.tight_layout(); fig.savefig(out / "fig_a3_hysteresis.png", dpi=200); plt.close(fig)

    # (b) phase diagram vs Ozbudak Fig 2c
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    Gs = np.logspace(-1, 3.3, 200)
    th = np.vstack([rn.thresholds_tmg(0.0, p)] + [rn.thresholds_tmg(G, p) for G in Gs])
    Gplot = np.r_[0.03, Gs]
    ax.fill_between(Gplot, th[:, 0], th[:, 1], color=GRID, alpha=0.9, lw=0, label="Model: two stable states")
    ax.plot(Gplot, th[:, 1], color=BLUE, lw=1.5); ax.plot(Gplot, th[:, 0], color=ORANGE, lw=1.5)
    d = FIG2C.astype(float); d.loc[d.glucose_uM == 0, "glucose_uM"] = 0.03
    ax.scatter(d.glucose_uM, d.switch_on_uM, s=46, color=BLUE, zorder=3, edgecolor="white", lw=1.5,
               label="Paper switch-ON (used to fit)")
    ax.scatter(d.glucose_uM, d.switch_off_uM, s=46, facecolor="white", edgecolor=ORANGE, lw=1.8, zorder=3,
               label="Paper switch-OFF (held out)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([0.03, 0.1, 1, 10, 100, 1000]); ax.set_xticklabels(["0", "0.1", "1", "10", "100", "1000"])
    ax.set_xlabel("Glucose outside the cell (µM)"); ax.set_ylabel("TMG outside the cell (µM)")
    ax.set_title("Model reproduces Ozbudak's phase diagram (Fig. 2c)", loc="left", fontsize=11.5)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "fig_a3_phase_diagram.png", dpi=200); plt.close(fig)

    # (c) rho estimates incl. Kuhlman
    rho_tab = pd.read_csv(out / "rho_estimates_rnaseq.csv")
    rows = [(f"{r.gene} RNA, vs {'same-batch WT' if 'same' in r.uninduced_reference else 'PRECISE-1K WT'}",
             r.rho, r.ci_lo, r.ci_hi) for r in rho_tab.itertuples()]
    rows += [("Fig. 2c thresholds (free fit)", res["rho_free_fit_to_fig2c"], np.nan, np.nan),
             ("Ozbudak 2004 (GFP protein)", 170, 140, 200),
             ("Kuhlman 2007, no cAMP (ΔlacY β-gal)", 258, np.nan, np.nan),
             ("Kuhlman 2007, 1 mM cAMP (ΔlacY β-gal)", 1650, np.nan, np.nan),
             ("Model, crp deleted (CRP-assisted loop lost)", rn.rho_eff(0.0, p, crp_knockout=True), np.nan, np.nan)]
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    for k, (lab, v, lo, hi) in enumerate(rows):
        c = INK if ("Ozbudak" in lab or "Fig" in lab) else AQUA if "Model" in lab else RED if "Kuhlman" in lab else BLUE
        if not np.isnan(lo):
            ax.plot([lo, hi], [k, k], color=c, lw=1.5)
        ax.scatter(v, k, color=c, s=40, zorder=3)
    ax.axvspan(1, 9, color=ORANGE, alpha=0.12, lw=0); ax.axvline(9, color=ORANGE, lw=1)
    ax.text(8.5, -0.4, "no switch\n(ρ < 9)", ha="right", va="bottom", color=ORANGE, fontsize=9)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    ax.set_xscale("log"); ax.set_xlim(3, 3000); ax.grid(axis="y", visible=False)
    ax.set_xlabel("Repression factor ρ (fully ON ÷ OFF)")
    ax.set_title("Every estimate of ρ is far above the switch limit of 9", loc="left", fontsize=11.5)
    fig.tight_layout(); fig.savefig(out / "fig_a3_rho.png", dpi=200); plt.close(fig)

    # (d) titration strains
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, rho in zip(axes, [170, 50, 5]):
        up, down = strains[rho]
        ax.loglog(Tg, up, color=BLUE, label="start OFF"); ax.loglog(Tg, down, color=ORANGE, ls=(0, (4, 2)), label="start ON")
        kind = "switch (memory)" if res[f"strain_rho{rho}"]["hysteretic"] else "dimmer (no memory)"
        ax.set_title(f"ρ = {rho}: {kind}", loc="left", fontsize=11); ax.set_xlabel("TMG (µM), 1 mM glucose")
    axes[0].set_ylabel("LacY per cell"); axes[0].legend(loc="upper left", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "fig_a3_titration.png", dpi=200); plt.close(fig)

    import matplotlib.image as mpimg
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))
    for ax, f in zip(axs.flat, ["fig_a3_hysteresis", "fig_a3_phase_diagram", "fig_a3_rho", "fig_a3_titration"]):
        ax.imshow(mpimg.imread(out / f"{f}.png")); ax.axis("off")
    fig.tight_layout(); fig.savefig(out / "hysteresis_check.png", dpi=150, facecolor="white"); plt.close(fig)


# ------------------------------------------------------------------ main
def main(net_csv, out="outputs/a3"):
    out = Path(out)
    p, fitted, outputs = load_inputs(out)
    res = {"version": "A3 v3 (network v3: A2 + RegulonDB 14.5.1 + Kuhlman 2007 + Ozbudak 2004)",
           "rho_free_fit_to_fig2c": fitted.get("rho_free_fit_to_fig2c")}

    # wild-type checks (must be unchanged from A3: the CRP-free basal only acts in a crp knockout)
    def rmse(col, k):
        pred = np.array([rn.thresholds_tmg(G, p)[k] for G in FIG2C.glucose_uM])
        return float(np.sqrt(np.mean(np.log10(pred / FIG2C[col].values) ** 2)))
    res["ozbudak_fig2c"] = {"switch_on_rmse_log10_fitted": rmse("switch_on_uM", 1),
                            "switch_off_rmse_log10_heldout": rmse("switch_off_uM", 0)}
    Ts = np.logspace(-0.3, 3.3, 90)
    sweeps = {}
    for G in [0.0, 1000.0]:
        up, down = sweep(p, G, Ts)
        off, on = switch_points(Ts, up, down, p, G)
        an_off, an_on = rn.thresholds_tmg(G, p)
        sweeps[G] = (up, down)
        res[f"ode_sweep_G{int(G)}"] = {"switch_off_ode_20h": float(off), "switch_on_ode_20h": float(on),
                                       "switch_off_analytic": float(an_off), "switch_on_analytic": float(an_on),
                                       "on_off_fold_expression": float(down.max() / up.min())}
    rhos = np.linspace(5, 15, 1001)
    res["cusp_rho"] = float(rhos[np.argmax([rn.switching_drives(r, p["n"]) is not None for r in rhos])])
    Tg = np.logspace(0, 4, 120)
    strains = {}
    for rho in [170, 50, 5]:
        up, down = sweep(p, 1000.0, Tg, rho=rho)
        strains[rho] = (up, down)
        res[f"strain_rho{rho}"] = {"hysteretic": bool(np.max(np.abs(np.log(down / up))) > 0.5),
                                   "hill_coefficient_up": float(hill_coefficient(Tg, up))}

    # crp knockout (Kuhlman: no cAMP response, ~1/240 basal; CRP-assisted looping lost)
    ko_up = np.array([rn.simulate(T, 0.0, p, induced=False, crp_knockout=True) for T in Ts])
    ko_th = rn.thresholds_tmg(0.0, p, crp_knockout=True)
    res["crp_knockout"] = {"rho_eff": rn.rho_eff(0.0, p, crp_knockout=True),
                           "crp_free_basal_fraction": p["crp_free_basal"],
                           "max_LacY_any_TMG_up_to_2000uM": float(ko_up.max()),
                           "wt_max_LacY": float(sweeps[0.0][1].max()),
                           "thresholds_tmg_uM": list(map(float, ko_th)) if ko_th else None,
                           "note": "rho stays > 9, so the loop is still bistable in principle, but the switch-on TMG "
                                   "moves far beyond any physiological dose: crp-deleted cells stay OFF (Traxler 2006)"}

    # ---------------- outputs
    circuit, context = edge_tables(Path(net_csv), outputs, Path(__file__).resolve().parents[1] / "data" / "curated" / "regulondb_crp_edges.json")
    res["edges"] = {"in_circuit_simulated": len(circuit), "context_not_simulated": len(context)}
    q = {k: (v if isinstance(v, str) or v is None or isinstance(v, bool) else float(v)) for k, v in p.items()}

    model = {
        "name": "crp_lac_circuit_A3_v3",
        "version": "v3 (24 Nov 2026)",
        "description": "Final parameterized CRP/lac circuit for A4. Lac core (Ozbudak 2004 form, fitted to his Fig 2c) "
                       "+ CRP-assisted LacI looping and CRP-free basal (Kuhlman 2007) + 21 A2-validated CRP outputs. "
                       "Edges carry their evidence (validated_by). The same model is in circuit_model.xml (SBML).",
        "changes_from_v2": [
            "crp_free_basal = 1/240 (Kuhlman Table 2): a crp knockout keeps a small lac output instead of exactly 0",
            "edges re-labelled with validated_by from network v3 (A2, RegulonDB 14.5.1 grade, Kuhlman intervention, Ozbudak)",
            "context_edges: RegulonDB regulators around the circuit (not simulated)",
            "SBML rewritten for COPASI + GillesPy2 (molecule counts, minutes), checked in both",
        ],
        "units": {"time": "min", "species": "molecules per cell", "TMG_ext": "uM", "glucose": "uM"},
        "inputs": {"TMG_ext": "extracellular TMG (uM)", "glucose": "extracellular glucose (uM)",
                   "crp_knockout": "bool: crp (or cyaA) deleted"},
        "species": {sid: {"description": doc,
                          "initial_uninduced": initial_counts(p, 0.0, False)[sid],
                          "initial_induced": initial_counts(p, 0.0, True)[sid]} for sid, doc in SPECIES_DOC.items()},
        "initial_conditions_note": "uninduced = Ozbudak's OFF history (grown without TMG), induced = ON history "
                                   "(pre-grown fully induced); values at 0 glucose. Use rn.initial_state(p, induced, G, crp_knockout) for others.",
        "reactions": [{"id": rid, "reactants": reac, "products": prod, "propensity": law}
                      for rid, reac, prod, law in REACTIONS],
        "rate_constants_by_condition": {
            "TMG15_G0": rate_constants(p, 15.0, 0.0),
            "TMG15_G1000": rate_constants(p, 15.0, 1000.0),
            "TMG15_G0_crp_knockout": rate_constants(p, 15.0, 0.0, crp_knockout=True),
        },
        "input_functions": {
            "k_rep_max": "(rho_eff - 1) * k_unrep",
            "k_tx_eff": "k_tx * c_eff;  c_eff = c_crp(glucose), or crp_free_basal if crp_knockout",
            "c_crp": "crp_min + (1 - crp_min) / (1 + glucose / crp_K)",
            "k_up_eff": "k_up * beta_T(TMG_ext) * beta_G(glucose)",
            "beta_T": "TMG_ext ** uptake_h   (power law, uM)",
            "beta_G": "excl_min + (1 - excl_min) / (1 + glucose / excl_K)",
            "d_y": "ln(2) / doubling_min",
            "rho_eff": "rho; crp knockout: 1 + (rho - 1) * (1 + L0) / (1 + L(crp_occ_ref)), "
                       "L(theta) = loop_L0 * (1 + (loop_omega - 1) * theta)  (Kuhlman CRP-assisted looping)",
        },
        "parameters": q,
        "parameter_sources": {
            "rho": "Ozbudak 2004 GFP (170); free fit to Fig 2c 173; PRECISE-1K lacY RNA 160",
            "n": "assumption (gives the cusp at rho = 9)",
            "k_up, crp_K, excl_K": "fitted to Ozbudak Fig 2c switch-on thresholds",
            "crp_min, excl_min, uptake_h": "Ozbudak Fig 3d-e",
            "loop_L0, loop_omega": "Kuhlman 2007 Fig 2a refit (paper Omega 10.3)",
            "crp_free_basal": "Kuhlman 2007 Table 2: 1 / f_cAMP = 1/240",
            "k_tx, d_m, k_tl, doubling_min, d_x, x0, k_unrep": "set (units / noise; A5 recalibrates)",
        },
        "derived": {"Y_max_LacY": rn.y_max(p), "basal_LacY": rn.y_max(p) / p["rho"],
                    "rho_by_condition": {"no_glucose": rn.rho_eff(0.0, p), "glucose_1mM": rn.rho_eff(1000.0, p),
                                         "crp_knockout": rn.rho_eff(0.0, p, crp_knockout=True)}},
        "validation": {k: res[k] for k in ["ozbudak_fig2c", "cusp_rho", "ode_sweep_G0", "ode_sweep_G1000",
                                           "strain_rho170", "strain_rho50", "strain_rho5", "crp_knockout"]},
        "edges": circuit,
        "context_edges": context,
        "crp_outputs": {
            "description": "A2-validated CRP targets: mRNA (TPM) = basal + amplitude * Hill(c), repressed targets 1 - Hill; "
                           "c = c_crp(glucose) (0 for crp knockout). Molecules = TPM * tpm_to_molecules.",
            "tpm_to_molecules": 0.0024, "targets": outputs},
        "files": {"sbml": "circuit_model.xml (TMG 15 uM, glucose 0, uninduced start; other conditions: "
                          "build_a3_v3.write_sbml(path, p, T, G, outputs, induced, crp_knockout))"},
    }
    model = json.loads(json.dumps(model, default=float))
    (out / "circuit_model.yaml").write_text(yaml.safe_dump(model, sort_keys=False, allow_unicode=True, width=120))

    ve = {"description": "Every edge of the A3 v3 circuit and where its validation comes from. "
                         "a2_validated = our own causal test on PRECISE-1K; validated_by lists all sources.",
          "sources": {"A2": "causal-mediation test on PRECISE-1K (a2_mediation_report.md)",
                      "Kuhlman2007_intervention": "graded cAMP / IPTG titration with crp-null control, PNAS 104:6043",
                      "RegulonDB_Confirmed_grade1": "RegulonDB 14.5.1: binding-site mutation + binding + expression",
                      "RegulonDB_<level>": "RegulonDB 14.5.1 confidence level",
                      "RegulonDB_grade<1-4>": "causal grade from RegulonDB evidence type (1 site mutation, 2 effect + binding, "
                                              "3 effect only, 4 binding only); only grades 1-2 with agreeing sign count",
                      "Ozbudak2004_measured": "measured in Ozbudak et al. 2004 Nature 427:737 (Fig 3)"},
          "edges": circuit}
    (out / "validated_edges.json").write_text(json.dumps(json.loads(json.dumps(ve, default=float)), indent=1,
                                                         allow_nan=False))

    sbml = out / "circuit_model.xml"
    write_sbml(sbml, p, 15.0, 0.0, outputs=outputs)
    res["sbml_check_TMG15_G0_uninduced"] = check_sbml(sbml, p, 15.0, 0.0)
    tmp = out / "_check_TMG40.xml"                                  # a run that actually switches ON
    write_sbml(tmp, p, 40.0, 0.0)
    res["sbml_check_TMG40_G0_uninduced"] = check_sbml(tmp, p, 40.0, 0.0)
    tmp.unlink()

    figures(p, out, res, sweeps, strains, Ts, Tg, (ko_up,))
    (out / "a3_checks_v3.json").write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main(*sys.argv[1:])
