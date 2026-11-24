"""
A6: Lamoureux 2023 Fig. 6 workflow applied to every implemented perturbation.

For each entry of config/precise1k_condition_map.yaml that has a matched real condition:

  real      control and perturbed samples -> PRECISE-1K activities (own A matrix, or projected like the AAT data)
            -> DiMA (perturbed - control), null from replicate pairs, BH FDR.
  simulated the calibrated circuit (A5) is run for the in-silico control and the perturbation (replicate populations).
            Each population becomes a pseudo-bulk sample: mean mRNA per cell of every circuit gene -> TPM, written into the
            real control condition's mean expression profile (all other genes unchanged); deleted genes set to 0.
            These samples are projected onto PRECISE-1K exactly like new data in Fig. 6 and compared by the same DiMA null.
  ceiling   the real perturbed values of the circuit genes only, written into the real control profile: the most any
            model of these genes could recover through the knowledge base.

Outputs: outputs/perturbation_dima_comparison.csv (one row per perturbation x iModulon)
         outputs/a6/perturbation_summary.csv, outputs/a6/simulated_samples_log_tpm.csv, outputs/a6/a6_results.json
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import reaction_network as rn  # noqa: E402
from a6_reproduce_aat import load_k12  # noqa: E402
from dima import P1K_DIR, P1K_REFS, KnowledgeBase, Panel, bh  # noqa: E402
from perturbation import load_condition_map, load_perturbations  # noqa: E402
from ssa_simulator import StochasticCircuit  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "outputs"
A6 = OUT / "a6"
MODEL = OUT / "calibrated_circuit_model.yaml"
LAC_GENES = {"lacZ": "b0344", "lacY": "b0343", "lacA": "b0342"}
DELETABLE = {"crp": "b3357", "cyaA": "b3806", "lacI": "b0345", "lacY": "b0343"}
KEY_IMODULONS = ["Crp-1", "Crp-2", "crp KO-1", "crp KO-2", "Galactose", "Maltose", "Ribose", "FlhDC-1", "FlhDC-2",
                 "Cra", "Thr/Ser-1", "Ferric Citrate", "Gluconate", "Fur-2"]


# ---------------------------------------------------------------------- knowledge base with every real condition
def build_kb(cmap):
    kb = KnowledgeBase.load()
    k12_projects = sorted({m[k].split(":")[0] for s in cmap["perturbations"] if (m := s.get("matched"))
                           and m["dataset"] == "k12_projected" for k in ("control", "perturbed")})
    lt, md = load_k12(kb, k12_projects)
    kb.add_samples(lt, md)
    kb.k12_log_tpm = lt
    kb.null_cdfs()                       # fixed now: simulated samples never enter the null
    return kb


def real_profile(kb, key):
    ids = kb.samples(key)
    src = [kb.log_tpm[i] if i in kb.log_tpm.columns else kb.k12_log_tpm[i] for i in ids]
    return pd.concat(src, axis=1).mean(axis=1), ids


# ---------------------------------------------------------------------- simulated pseudo-bulk samples
class PseudoBulk:
    """Simulated population -> log2(TPM+1) for the circuit genes, in PRECISE-1K units."""

    def __init__(self, sc: StochasticCircuit, kb: KnowledgeBase):
        self.sc, self.kb = sc, kb
        gt = kb.gene_table
        self.name2b = {n: b for b, n in gt.gene_name.items()}
        self.targets = {g: self.name2b[g] for g in sc.outputs}
        # lac genes: anchor molecules -> TPM once, at the knowledge-base reference (WT, M9 glucose 2 g/L, no inducer),
        # using the analytic OFF-state mRNA mean (the sampled mean of ~0.01 molecules per cell is too noisy to anchor on)
        p = sc.params
        G = 2.0 * 1e6 / 180.16
        m_ref = rn.initial_state(p, False, G, False)[1]
        ref_tpm = 2 ** kb.log_tpm.loc[list(LAC_GENES.values()), P1K_REFS].mean(axis=1) - 1
        self.lac_scale = {b: float(ref_tpm[b] / m_ref) for b in LAC_GENES.values()}
        self.tpm_to_molecules = sc.tpm_to_molecules

    def sample(self, pert, seed, n_cells):
        tr = self.sc.run(pert, seed=seed, n_cells=n_cells)
        snap = tr.at(pert.t_end_min)
        tpm = {b: snap[f"mRNA_{g}_count"].mean() / self.tpm_to_molecules for g, b in self.targets.items()}
        for b in LAC_GENES.values():
            tpm[b] = snap.M_count.mean() * self.lac_scale[b]
        x = np.log2(pd.Series(tpm) + 1)
        for g in pert.deleted_genes:
            if g in DELETABLE:
                x[DELETABLE[g]] = 0.0               # gene removed from the genome: no reads
        return x, {"Y_mean": float(snap.Y_count.mean()), "frac_ON": float((snap.Y_count > tr.meta["on_threshold_Y"]).mean())}


def hybrid(base: pd.Series, circuit: pd.Series):
    x = base.copy()
    x.loc[circuit.index] = circuit.values
    return x


# ---------------------------------------------------------------------- DiMA on given activity vectors
def dima_vectors(kb, a1: pd.DataFrame, a2: pd.DataFrame):
    null = kb.null_cdfs()
    d = a2.mean(axis=1) - a1.mean(axis=1)
    p = pd.Series({k: 1 - null[k].cdf(abs(d[k])) for k in d.index})
    return d, p, bh(p)


def sig(d, q, s):
    return (q < s["fdr"]) & (d.abs() > s["threshold"])


def thresholds():
    return json.load(gzip.open(P1K_DIR / "p1k_precise1k.json.gz"))["thresholds"]


def membership(kb):
    th = pd.Series(thresholds())
    return {im: set(kb.M.index[kb.M[im].abs() > th[im]]) for im in kb.M.columns}


EXCLUDE_SYSTEMS = ("Genetic Alterations", "ALE Effects", "Technical")   # strain / evolution / artefact iModulons


def build_panel(kb, circ_b):
    P = Panel(kb, sorted(circ_b), thresholds(), min_fidelity=0.5)
    sysc = kb.imodulon_table["system_category"].astype(str)
    P.readable = [im for im in P.readable if not sysc.get(im, "").startswith(EXCLUDE_SYSTEMS)]
    return P


def main(n_cells=2000, reps=None, model=MODEL, a6_dir=A6, main_csv=None, only=None):
    """model: circuit model to simulate; a6_dir: where per-run tables go; main_csv: path of the main comparison table;
    only: optional list of perturbation names to run."""
    A6 = Path(a6_dir); A6.mkdir(parents=True, exist_ok=True)
    main_csv = Path(main_csv) if main_csv else OUT / "perturbation_dima_comparison.csv"
    cmap = load_condition_map()
    S = cmap["dima_settings"]
    reps = reps or cmap["defaults"].get("replicates", 3)
    kb = build_kb(cmap)
    sc = StochasticCircuit(model)
    pb = PseudoBulk(sc, kb)
    mem = membership(kb)
    circ_b = set(pb.targets.values()) | set(LAC_GENES.values()) | set(DELETABLE.values())
    touched = {im for im, g in mem.items() if g & circ_b}
    panel = build_panel(kb, circ_b)
    prow, grow = [], []

    rows, summary, sim_cols, extra = [], [], {}, {}
    for pert, ctrl, m in load_perturbations(with_controls=True):
        if not m or (only and pert.name not in only):
            continue
        base, ctrl_ids = real_profile(kb, m["control"])
        pert_prof, pert_ids = real_profile(kb, m["perturbed"])
        # simulated replicate populations
        X_ctrl, X_pert, info = {}, {}, {}
        for r in range(1, reps + 1):
            xc, ic = pb.sample(ctrl, seed=100 + r, n_cells=n_cells)
            xp, ip = pb.sample(pert, seed=200 + r, n_cells=n_cells)
            X_ctrl[f"sim__{pert.name}__control__{r}"] = hybrid(base, xc)
            X_pert[f"sim__{pert.name}__perturbed__{r}"] = hybrid(base, xp)
            info[r] = {"control": ic, "perturbed": ip}
        Xc, Xp = pd.DataFrame(X_ctrl), pd.DataFrame(X_pert)
        Ac, Ap = kb.infer_activities(kb.center(Xc)), kb.infer_activities(kb.center(Xp))
        sim_d, sim_p, sim_q = dima_vectors(kb, Ac, Ap)
        # real
        real = kb.dima(m["control"], m["perturbed"], **S["fig6"])
        real_d, real_q = real.difference, real.qvalue
        a_ctrl_real = kb.A[ctrl_ids].mean(axis=1)
        # ceiling: real circuit-gene changes only
        cg = sorted(circ_b & set(base.index))
        Xceil = pd.DataFrame({"ceil_c": base, "ceil_p": hybrid(base, pert_prof.loc[cg])})
        Aceil = kb.infer_activities(kb.center(Xceil))
        ceil_d = Aceil["ceil_p"] - Aceil["ceil_c"]
        ceil_p = pd.Series({k: 1 - kb.null_cdfs()[k].cdf(abs(ceil_d[k])) for k in ceil_d.index}); ceil_q = bh(ceil_p)

        for im in kb.A.index:
            rows.append({"perturbation": pert.name, "kind": pert.kind, "imodulon": im,
                         "matched_control": m["control"], "matched_perturbed": m["perturbed"], "dataset": m["dataset"],
                         "held_out": m["held_out"],
                         "real_activity_control": a_ctrl_real[im], "real_activity_perturbed": real.loc[im, "activity_2"],
                         "real_difference": real_d[im], "real_qvalue": real_q[im],
                         "real_dima_fig6": bool(sig(real_d, real_q, S["fig6"])[im]),
                         "real_dima_default": bool(sig(real_d, real_q, S["default"])[im]),
                         "sim_activity_control": a_ctrl_real[im],
                         "sim_activity_perturbed": a_ctrl_real[im] + sim_d[im],
                         "sim_difference": sim_d[im], "sim_qvalue": sim_q[im],
                         "sim_dima_fig6": bool(sig(sim_d, sim_q, S["fig6"])[im]),
                         "sim_dima_default": bool(sig(sim_d, sim_q, S["default"])[im]),
                         "ceiling_difference": ceil_d[im], "ceiling_dima_default": bool(sig(ceil_d, ceil_q, S["default"])[im]),
                         "contains_circuit_gene": im in touched})
        # circuit-gene panel readout: same samples, activities read only from the circuit's genes
        real_c = pd.concat([kb.log_tpm[i] if i in kb.log_tpm.columns else kb.k12_log_tpm[i] for i in ctrl_ids], axis=1)
        real_p = pd.concat([kb.log_tpm[i] if i in kb.log_tpm.columns else kb.k12_log_tpm[i] for i in pert_ids], axis=1)
        Pr_c, Pr_p = panel.activities(kb.center(real_c)), panel.activities(kb.center(real_p))
        Ps_c, Ps_p = panel.activities(kb.center(Xc)), panel.activities(kb.center(Xp))
        pr_d, _, pr_q = panel.dima(Pr_c, Pr_p)
        ps_d, _, ps_q = panel.dima(Ps_c, Ps_p)
        for im in panel.readable:
            prow.append({"perturbation": pert.name, "imodulon": im, "panel_genes": ";".join(kb.gene_table.gene_name[panel.members[im]]),
                         "panel_fidelity": panel.fidelity[im], "held_out": m["held_out"],
                         "panel_real_control": Pr_c.loc[im].mean(), "panel_real_perturbed": Pr_p.loc[im].mean(),
                         "panel_real_difference": pr_d[im], "panel_real_qvalue": pr_q[im],
                         "panel_sim_difference": ps_d[im], "panel_sim_qvalue": ps_q[im],
                         "panel_real_dima": bool(sig(pr_d, pr_q, S["default"])[im]),
                         "panel_sim_dima": bool(sig(ps_d, ps_q, S["default"])[im])})
        # gene level: log2 fold change of every circuit gene (what the simulator outputs directly)
        for b in cg:
            grow.append({"perturbation": pert.name, "gene": kb.gene_table.gene_name[b], "b_number": b,
                         "held_out": m["held_out"],
                         "real_log2fc": float(real_p.loc[b].mean() - real_c.loc[b].mean()),
                         "real_log2fc_sd": float(np.sqrt(real_p.loc[b].var(ddof=1) / len(pert_ids) + real_c.loc[b].var(ddof=1) / len(ctrl_ids))),
                         "sim_log2fc": float(Xp.loc[b].mean() - Xc.loc[b].mean()),
                         "real_control_log_tpm": float(real_c.loc[b].mean()), "sim_control_log_tpm": float(Xc.loc[b].mean()),
                         "is_deleted_gene": any(DELETABLE.get(g) == b for g in pert.deleted_genes)})
        sim_cols.update({**X_ctrl, **X_pert})
        extra[pert.name] = {"populations": info, "real_samples": {"control": ctrl_ids, "perturbed": pert_ids}}

    df = pd.DataFrame(rows)
    df.to_csv(main_csv, index=False)
    pd.DataFrame(sim_cols).loc[sorted(circ_b & set(pd.DataFrame(sim_cols).index))].to_csv(A6 / "simulated_samples_log_tpm.csv")

    for name, d in df.groupby("perturbation", sort=False):
        r = {"perturbation": name, "matched": f"{d.matched_control.iloc[0]} -> {d.matched_perturbed.iloc[0]}",
             "held_out": bool(d.held_out.iloc[0])}
        for lvl in ("fig6", "default"):
            real_s, sim_s = d[f"real_dima_{lvl}"], d[f"sim_dima_{lvl}"]
            agree = real_s & sim_s & (np.sign(d.real_difference) == np.sign(d.sim_difference))
            r[f"real_dimas_{lvl}"] = int(real_s.sum())
            r[f"sim_dimas_{lvl}"] = int(sim_s.sum())
            r[f"recovered_{lvl}"] = int(agree.sum())
            r[f"recovered_names_{lvl}"] = ";".join(d.imodulon[agree])
        ceil_hit = d.real_dima_default & d.ceiling_dima_default & (np.sign(d.real_difference) == np.sign(d.ceiling_difference))
        r["ceiling_recovered_default"] = int(ceil_hit.sum())
        r["r_all_imodulons"] = float(np.corrcoef(d.sim_difference, d.real_difference)[0, 1])
        t = d[d.contains_circuit_gene]
        r["r_circuit_imodulons"] = float(np.corrcoef(t.sim_difference, t.real_difference)[0, 1])
        r["r_sim_vs_ceiling"] = float(np.corrcoef(d.sim_difference, d.ceiling_difference)[0, 1])
        r["ceiling_r_vs_real"] = float(np.corrcoef(d.ceiling_difference, d.real_difference)[0, 1])
        k = d[d.contains_circuit_gene & (d.real_difference.abs() > 2)]
        r["sign_agreement_circuit_imodulons"] = float((np.sign(k.sim_difference) == np.sign(k.real_difference)).mean()) if len(k) else None
        summary.append(r)
    Gdf = pd.DataFrame(grow)
    Gdf.to_csv(A6 / "gene_level_comparison.csv", index=False)
    P = pd.DataFrame(prow)
    P.to_csv(A6 / "panel_dima_comparison.csv", index=False)
    for r in summary:
        d = P[P.perturbation == r["perturbation"]]
        hit = d.panel_real_dima & d.panel_sim_dima & (np.sign(d.panel_real_difference) == np.sign(d.panel_sim_difference))
        r["panel_real_dimas"] = int(d.panel_real_dima.sum())
        r["panel_sim_dimas"] = int(d.panel_sim_dima.sum())
        r["panel_recovered"] = int(hit.sum())
        r["panel_recovered_names"] = ";".join(d.imodulon[hit])
        r["panel_false_positive_names"] = ";".join(d.imodulon[d.panel_sim_dima & ~d.panel_real_dima])
        g = Gdf[(Gdf.perturbation == r["perturbation"]) & ~Gdf.is_deleted_gene]
        moved = g[g.real_log2fc.abs() > 1]
        r["gene_r"] = float(np.corrcoef(g.sim_log2fc, g.real_log2fc)[0, 1])
        r["gene_sign_agreement_moved"] = float((np.sign(moved.sim_log2fc) == np.sign(moved.real_log2fc)).mean()) if len(moved) else None
        r["gene_n_moved"] = int(len(moved))
        strict = sig(d.panel_real_difference, d.panel_real_qvalue, S["fig6"]), sig(d.panel_sim_difference, d.panel_sim_qvalue, S["fig6"])
        hit_s = strict[0] & strict[1] & (np.sign(d.panel_real_difference) == np.sign(d.panel_sim_difference))
        r["panel_fig6_real_dimas"], r["panel_fig6_sim_dimas"], r["panel_fig6_recovered"] = int(strict[0].sum()), int(strict[1].sum()), int(hit_s.sum())
        multi = d.panel_genes.str.count(";") >= 1
        hit_m = hit & multi
        r["panel_multigene_real_dimas"], r["panel_multigene_sim_dimas"], r["panel_multigene_recovered"] = \
            int((d.panel_real_dima & multi).sum()), int((d.panel_sim_dima & multi).sum()), int(hit_m.sum())
        r["panel_r"] = float(np.corrcoef(d.panel_sim_difference, d.panel_real_difference)[0, 1])
        big = d[d.panel_real_difference.abs() > 2]
        r["panel_sign_agreement"] = float((np.sign(big.panel_sim_difference) == np.sign(big.panel_real_difference)).mean()) if len(big) else None
    pd.DataFrame(summary).to_csv(A6 / "perturbation_summary.csv", index=False)
    res = {"n_cells_per_population": n_cells, "replicates": reps, "circuit_genes_b": sorted(circ_b),
           "imodulons_containing_circuit_genes": sorted(touched), "lac_tpm_anchor": pb.lac_scale,
           "panel": {"readable": panel.readable, "fidelity": panel.fidelity.round(3).to_dict(),
                     "members": {k: list(kb.gene_table.gene_name[v]) for k, v in panel.members.items()}},
           "summary": summary, "populations": extra}
    (A6 / "a6_results.json").write_text(json.dumps(res, indent=2, default=float))
    return df, pd.DataFrame(summary)


if __name__ == "__main__":
    df, s = main()
    pd.set_option("display.width", 250)
    print(s[["perturbation", "held_out", "real_dimas_fig6", "sim_dimas_fig6", "recovered_fig6", "real_dimas_default",
             "sim_dimas_default", "recovered_default", "ceiling_recovered_default", "r_all_imodulons", "r_circuit_imodulons",
             "r_sim_vs_ceiling", "sign_agreement_circuit_imodulons"]].round(2).to_string())
    print(s[["perturbation", "gene_r", "gene_n_moved", "gene_sign_agreement_moved", "panel_real_dimas", "panel_sim_dimas", "panel_recovered", "panel_r", "panel_sign_agreement"]].round(2).to_string())
    for _, r in s.iterrows():
        print(r.perturbation, "| panel recovered:", r.panel_recovered_names, "| sim-only:", r.panel_false_positive_names)
