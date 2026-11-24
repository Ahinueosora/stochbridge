"""
A5 validation: does the calibrated simulator produce PETRI-seq-like cells?

  1. Simulate exponential (EZ Rich, glucose 11.1 mM) and stationary (glucose exhausted) cells with the A4 engine on
     outputs/calibrated_circuit_model.yaml (and, for contrast, the uncalibrated A3 v3 model), same cell numbers as the data.
  2. Pass molecules through the observation model: UMI ~ Poisson(q_g * kappa_i * molecules), kappa_i from the phase's
     real total-mRNA-UMI distribution; q_g fitted on EXPONENTIAL means only.
  3. Compare with the real cells:
       exponential  - shape only (Fano, zero fraction, histogram): means are matched by q_g
       stationary   - everything is a prediction (means through CRP activity 0.2 -> 1 and capture; shapes)
       phase split  - the same 17-operon PCA as on the real data (AUC, PC1 histograms)
  4. A4 memory check again with the calibrated lac bursts (maintenance dose, both histories).

Usage: python src/a5_validate.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from fit_noise_params import G_EXP, G_STAT, capture_sampler, load_summary  # noqa: E402
from load_scrnaseq import MODEL_OPERONS, OPERON_TO_MODEL_OUTPUT, auc_and_balacc, counts_pca  # noqa: E402
from perturbation import InducerLevelChange, Perturbation  # noqa: E402
from ssa_simulator import StochasticCircuit  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "outputs" / "a5"
OPS = list(MODEL_OPERONS)                     # lacZYA + 16 CRP-output operons


def simulate_phase(sc, phase, n, seed, overrides=None):
    G = G_EXP if phase == "exp" else G_STAT
    pt = Perturbation(name=f"petri_{phase}", TMG_uM=0.0, glucose_uM=G, n_cells=n, t_end_min=600, dt_min=600, seed=seed,
                      param_overrides=overrides or {})
    tr = sc.run(pt)
    snap = tr.at(600)
    mol = {"lacZYA": snap.M_count.values}
    for op, g in OPERON_TO_MODEL_OUTPUT.items():
        mol[op] = snap[f"mRNA_{g}_count"].values
    return pd.DataFrame(mol)[OPS], tr


def observe(mol, tot, q, rng):
    kappa = tot / 222.0
    return pd.DataFrame({op: rng.poisson(q[op] * kappa * mol[op].values) for op in OPS})


def stats(Y):
    out = {}
    for op in OPS:
        y = Y[op].values
        m = y.mean()
        out[op] = {"mean": float(m), "fano": float(y.var(ddof=1) / m) if m > 0 else None,
                   "frac0": float((y == 0).mean()), "hist": [int((y == k).sum()) for k in range(4)] + [int((y >= 4).sum())]}
    return out


def medium_scales(sc, S, q0):
    """Expression levels in EZ Rich differ from the PRECISE-1K fits: scale each gene's transcription so that, at the
    common capture q0, the simulated exponential mean equals the observed one (fitted on exponential cells only).
    Lac core: same for basal lac transcription (k_tx) - at TMG 0 cells cannot switch, so only the OFF level is affected."""
    pt = Perturbation(TMG_uM=0.0, glucose_uM=G_EXP)
    cond = sc.condition(pt)
    p = sc.effective_params(pt)
    scale = {}
    for op, g in OPERON_TO_MODEL_OUTPUT.items():
        m_model = cond["out_rates"][f"k_{g}"] / p["d_m"]                    # mean molecules per cell
        m_data = S["operons"][op]["exp"][0] / q0                             # kappa averages to 1 in exp
        scale[g] = m_data / m_model if m_model > 0 else 1.0
    m_lac = p["k_tx"] * cond["rates"]["k_tx_eff"] / p["k_tx"] / p["d_m"] / p["rho"]
    lac_scale = (S["operons"]["lacZYA"]["exp"][0] / q0) / m_lac
    return scale, lac_scale


def run_model(yaml_path, S, label, seed=11, q0=0.05):
    rng = np.random.default_rng(seed)
    sc = StochasticCircuit(yaml_path, mode="fast", include_outputs=True)
    n = {"exp": S["cells"]["n_exp"], "stat": S["cells"]["n_stat"]}
    scale, lac_scale = medium_scales(sc, S, q0)
    ov = {"output_scale": scale, "k_tx": sc.params["k_tx"] * lac_scale}
    mol, tot, trs = {}, {}, {}
    for k, ph in enumerate(("exp", "stat")):
        mol[ph], trs[ph] = simulate_phase(sc, ph, n[ph], seed + k, ov)
        tot[ph] = capture_sampler(S, ph, n[ph], rng)
    q = {op: q0 for op in OPS}                       # one capture per molecule for every gene
    Y = {ph: observe(mol[ph], tot[ph], q, rng) for ph in ("exp", "stat")}
    st = {ph: stats(Y[ph]) for ph in Y}
    # phase split on simulated cells, same pipeline as the real data
    C = np.vstack([Y["exp"].values, Y["stat"].values]).astype(float)
    T = np.r_[tot["exp"], tot["stat"]]
    is_exp = np.r_[np.ones(n["exp"], bool), np.zeros(n["stat"], bool)]
    s1, w = counts_pca(C, T, gm=S["cells"]["geomean_mrna_umi"])
    if w[OPS.index("fecE")] < 0:
        s1, w = -s1, -w
    auc, bal = auc_and_balacc(s1, is_exp)
    lo, hi = S["phase_split"]["model17_operons_pc1"]["pc1_hist_24bins"]["range"]
    hist = {ph: np.histogram(s1[m], bins=24, range=(min(lo, s1.min()), max(hi, s1.max())))[0].tolist()
            for ph, m in (("exp", is_exp), ("stat", ~is_exp))}
    return {"label": label, "yaml": str(Path(yaml_path).name), "q0": q0, "medium_scale": scale, "lac_basal_scale": lac_scale,
            "stats": st,
            "pca17": {"auc": auc, "balanced_accuracy": bal, "loadings": dict(zip(OPS, map(float, w))),
                      "pc1_scores_exp": s1[is_exp].tolist(), "pc1_scores_stat": s1[~is_exp].tolist(), "pc1_hist": hist},
            "molecules_mean": {ph: mol[ph].mean().to_dict() for ph in mol},
            "crp_outputs_method": trs["exp"].meta.get("crp_outputs_method")}


def compare(S, M):
    """Agreement metrics between simulated and real per-gene summaries."""
    rows = []
    for ph in ("exp", "stat"):
        for op in OPS:
            r = S["operons"][op][ph]
            s = M["stats"][ph][op]
            n = sum(r[2])
            rf = r[1] / r[0] if r[0] > 0 else None
            rows.append({"phase": ph, "operon": op, "real_mean": r[0], "sim_mean": s["mean"], "real_frac0": r[2][0] / n,
                         "sim_frac0": s["frac0"], "real_fano": rf, "sim_fano": s["fano"],
                         "real_hist": r[2], "sim_hist": s["hist"]})
    df = pd.DataFrame(rows)
    out = {}
    for ph in ("exp", "stat"):
        d = df[(df.phase == ph) & (df.real_mean > 0) & (df.sim_mean > 0)]
        out[ph] = {"log10_mean_r": float(np.corrcoef(np.log10(d.real_mean), np.log10(d.sim_mean))[0, 1]),
                   "log10_mean_rmse": float(np.sqrt(np.mean((np.log10(d.real_mean) - np.log10(d.sim_mean)) ** 2))),
                   "median_abs_log2_mean_ratio": float(np.median(np.abs(np.log2(d.sim_mean / d.real_mean))))}
        e = d[(d.real_mean >= 0.05)]
        if len(e):
            out[ph]["fano_rmse_expressed"] = float(np.sqrt(np.mean((e.real_fano - e.sim_fano) ** 2)))
            out[ph]["fano_pairs_expressed"] = {o: [float(a), float(b)] for o, a, b in zip(e.operon, e.real_fano, e.sim_fano)}
        # histogram chi2 over all genes (bins 0,1,2,3,4+ merged where expected < 5)
        chi, dof = 0.0, 0
        for _, rw in d.iterrows():
            ro, so = np.array(rw.real_hist, float), np.array(rw.sim_hist, float)
            so = so / so.sum() * ro.sum()
            keep = so >= 5
            if keep.sum() >= 2:
                ro2 = np.r_[ro[keep], ro[~keep].sum()] if (~keep).any() else ro[keep]
                so2 = np.r_[so[keep], so[~keep].sum()] if (~keep).any() else so[keep]
                ok = so2 > 0
                chi += float(np.sum((ro2[ok] - so2[ok]) ** 2 / so2[ok])); dof += int(ok.sum() - 1)
        out[ph]["hist_chi2"] = chi; out[ph]["hist_dof"] = dof
    return out, df


def memory_check(yaml_path, n=1000, seed=31):
    sc = StochasticCircuit(yaml_path, mode="fast", include_outputs=False)
    res = {}
    for hist in ("uninduced", "induced"):
        tr = sc.run(InducerLevelChange(11.0, initial_state=hist, n_cells=n, t_end_min=1800, dt_min=10, seed=seed))
        W = tr.wide("Y").values; thr = tr.meta["on_threshold_Y"]
        hit = (W > thr) if hist == "uninduced" else (W < thr)
        res[hist] = {"fraction_flipped_30h": float(hit.any(1).mean()), "fraction_ON_20h": float((tr.at(1200).Y_count > thr).mean()),
                     "basal_LacY_mean_t0_2h": float(W[:, :13].mean()), "basal_LacY_fano_2h": float(W[:, 12].var() / max(W[:, 12].mean(), 1e-9))}
    sweep = []
    for T in (5.85, 8.33, 11.85, 16.87, 24.02, 34.2, 48.69):
        tr = sc.run(InducerLevelChange(T, n_cells=n, seed=int(T * 10)))
        sweep.append({"TMG_uM": T, "fraction_ON_20h_from_OFF": float((tr.at(1200).Y_count > tr.meta["on_threshold_Y"]).mean())})
    res["sweep_from_OFF"] = sweep
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    S = load_summary()
    real = S["phase_split"]["model17_operons_pc1"]
    R = {"real": {"pca17_auc": real["auc"], "pca17_balanced_accuracy": real["balanced_accuracy"],
                  "all_operon_pc1_auc": S["phase_split"]["all_operons_pc1"]["auc"]}}
    M = {"calibrated": run_model(ROOT / "outputs" / "calibrated_circuit_model.yaml", S, "A5 calibrated"),
         "uncalibrated": run_model(ROOT / "outputs" / "circuit_model.yaml", S, "A3 v3 (uncalibrated)")}
    tables = []
    for k, m in M.items():
        R[k] = {kk: v for kk, v in m.items() if kk not in ("stats",)}
        R[k]["pca17"] = {kk: v for kk, v in m["pca17"].items() if not kk.startswith("pc1_scores")}
        R[k]["agreement"], df = compare(S, m)
        df["model"] = k; tables.append(df)
    pd.concat(tables).to_csv(OUT / "per_gene_comparison.csv", index=False)
    np.savez_compressed(OUT / "pc1_scores.npz", **{f"{k}_{ph}": np.array(M[k]["pca17"][f"pc1_scores_{ph}"]) for k in M for ph in ("exp", "stat")})
    R["memory_check"] = {"calibrated": memory_check(ROOT / "outputs" / "calibrated_circuit_model.yaml"),
                         "uncalibrated_A4": memory_check(ROOT / "outputs" / "circuit_model.yaml")}
    (OUT / "validation.json").write_text(json.dumps(R, indent=2, default=float))
    return R


if __name__ == "__main__":
    R = main()
    for k in ("calibrated", "uncalibrated"):
        print(k, "PCA17", R[k]["pca17"]["auc"], R[k]["pca17"]["balanced_accuracy"], "agreement", json.dumps(R[k]["agreement"], default=float)[:900])
    print("real", R["real"])
    print("memory", json.dumps(R["memory_check"], default=float))
