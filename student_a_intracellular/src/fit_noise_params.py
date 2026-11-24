"""
A5 - noise_calibration: fit the noise parameters of the A3 v3 circuit against PETRI-seq (GSE141018, exp 2.01) and the
single-molecule lac literature, and write outputs/calibrated_circuit_model.yaml.

Noise layers (what each source can and cannot constrain):
  1. Technical (observation) noise   PETRI-seq   per-cell capture kappa_i (empirical, per phase); capture per molecule q0
  2. Gene-level extrinsic noise     PETRI-seq   eta2: per-cell, per-gene CV2 of the transcription rate (slope of
                                                  Fano-1 vs mean, after removing capture variation)
  3. Transcriptional bursting       PETRI-seq   b: intercept of Fano-1 vs mean -> bounded only (5% capture hides it)
  4. Shared CRP-activity noise      PETRI-seq   sigma_c: upper bound from the near-zero correlations between CRP targets
  5. Lac protein burst size         Yu et al. 2006 (repressed lac promoter): 4.2 +/- 0.5 proteins per mRNA; PETRI-seq
                                                  cannot see it (mRNA only; lac is OFF: 0.013 UMIs per cell)
Method: weighted method of moments for 2-3, then an ABC check (simulate the observation model on a parameter grid,
accept parameter sets whose simulated Fano / zero-fraction / correlation summaries match the data).

Usage: python src/fit_noise_params.py [outputs dir]
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import reaction_network as rn  # noqa: E402
from load_scrnaseq import OPERON_TO_MODEL_OUTPUT  # noqa: E402

ROOT = HERE.parent
SUMMARY = ROOT / "data" / "scrnaseq" / "GSE141018" / "petri_2_01_summary.json"
G_EXP, G_STAT = 11100.0, 0.0      # EZ Rich glucose 0.2% = 11.1 mM; stationary: glucose exhausted
Q0, Q0_RANGE = 0.05, (0.025, 0.10)  # PETRI-seq capture per mRNA molecule (paper: ~2.5-10%)
YU2006 = {"burst_size": 4.2, "burst_size_se": 0.5, "bursts_per_cycle": 1.2, "cell_cycle_min": 55,
          "mRNA_lifetime_min": 1.5, "copies_per_cell": 5.0}
CHOI2008 = {"switch_threshold_LacY": 375, "uninduced_LacY_range": [0, 10]}


# ------------------------------------------------------------------ data
def load_summary(path=SUMMARY):
    return json.loads(Path(path).read_text())


def capture_sampler(S, phase, n, rng):
    """Draw per-cell total mRNA UMIs from the phase's empirical distribution (log2 half-bins, uniform within bin)."""
    h = np.array(S["mrna_umi"][phase]["hist_log2_halfbins"], float)
    b = rng.choice(len(h), n, p=h[:, 1] / h[:, 1].sum())
    return 2 ** (h[b, 0] + rng.uniform(0, 0.5, n))


def gene_moments(S, phase, min_mean=0.004):
    rows = []
    for g, r in S["operons"].items():
        if g.startswith("ctrl") or g.startswith("_") or r[phase][2] is None:
            continue
        mu, var = r[phase][0], r[phase][1]
        if mu >= min_mean:
            rows.append((g, mu, var / mu))
    return rows


# ------------------------------------------------------------------ 1-3: method of moments
def fit_fano_line(S, phase, cv2_kappa):
    """Fano-1 = A + B mu.  Observation model y ~ Poisson(q kappa m), m with CV2 (1+b)/M + eta2:
       B = cv2_kappa + eta2 (1 + cv2_kappa);  A = q0 (1 + b)(1 + cv2_kappa)."""
    n = S["cells"]["n_exp" if phase == "exp" else "n_stat"]
    rows = gene_moments(S, phase)
    mu = np.array([r[1] for r in rows]); F = np.array([r[2] for r in rows])
    se = np.sqrt((1 / mu + 2 * F ** 2) / n)
    W = 1 / se ** 2
    A = np.vstack([np.ones_like(mu), mu]).T
    beta = np.linalg.solve((A * W[:, None]).T @ A, (A * W[:, None]).T @ (F - 1))
    cov = np.linalg.inv((A * W[:, None]).T @ A)
    a, b = beta; sa, sb = np.sqrt(np.diag(cov))
    eta2 = (b - cv2_kappa) / (1 + cv2_kappa)
    out = {"phase": phase, "n_genes": len(rows), "intercept_A": a, "intercept_se": sa, "slope_B": b, "slope_se": sb,
           "capture_cv2": cv2_kappa, "eta2": eta2, "eta2_se": sb / (1 + cv2_kappa),
           "burst_b_by_q0": {str(q): {"b": a / (q * (1 + cv2_kappa)) - 1,
                                      "b_upper95": (a + 1.96 * sa) / (q * (1 + cv2_kappa)) - 1} for q in (Q0_RANGE[0], Q0, Q0_RANGE[1])},
           "genes": {g: {"mean": m_, "fano": f_} for g, m_, f_ in rows}}
    return out


# ------------------------------------------------------------------ 4: shared CRP noise bound
def elasticity(p, outputs, gene, G):
    """d ln(TPM_g) / d ln c at the condition's CRP activity."""
    c0 = rn.c_crp(G, p)
    def tpm(c):
        o = outputs[gene]
        a = c ** o["hill_n"] / (o["K_crp"] ** o["hill_n"] + c ** o["hill_n"])
        f = a if o["direction"] == "activated" else 1 - a
        return o["basal_TPM"] + o["amplitude_TPM"] * f
    h = 1e-3
    return (np.log(tpm(c0 * (1 + h))) - np.log(tpm(c0 * (1 - h)))) / (np.log(1 + h) - np.log(1 - h))


def crp_noise_bound(S, p, outputs):
    """Shared CRP-activity noise (lognormal sigma_c) makes CRP targets co-vary: for capture-scaled counts,
    cov(a,b) ~ mu_a mu_b e_a e_b sigma_c^2. Compare with the observed correlations (exp) and their SE = 1/sqrt(n)."""
    genes = S["capture_scaled_correlations"]["genes"]
    R = np.array(S["capture_scaled_correlations"]["exp"])
    n = S["cells"]["n_exp"]
    pairs = [("rbs", "galK"), ("rbs", "hflC"), ("galK", "hflC"), ("gyrA", "fecE"), ("rbs", "gyrA"), ("galK", "fecE"),
             ("hflC", "gyrA"), ("rbs", "fecE"), ("galK", "gyrA"), ("hflC", "fecE")]
    m = OPERON_TO_MODEL_OUTPUT
    rows = []
    for a, b in pairs:
        ea, eb = elasticity(p, outputs, m[a], G_EXP), elasticity(p, outputs, m[b], G_EXP)
        ma, mb = S["operons"][a]["exp"][0], S["operons"][b]["exp"][0]
        va, vb = S["operons"][a]["exp"][1], S["operons"][b]["exp"][1]
        k = ma * mb * ea * eb / np.sqrt(va * vb)            # predicted r per unit sigma_c^2
        rows.append({"pair": f"{a}-{b}", "elasticities": [ea, eb], "r_obs": float(R[genes.index(a), genes.index(b)]),
                     "r_per_sigma2": float(k)})
    k = np.array([r["r_per_sigma2"] for r in rows]); r = np.array([r["r_obs"] for r in rows])
    w = np.ones_like(k) * n                               # var(r) ~ 1/n
    s2 = float(np.sum(w * k * r) / np.sum(w * k * k))
    s2_se = float(1 / np.sqrt(np.sum(w * k * k)))
    return {"pairs": rows, "sigma_c2_hat": s2, "sigma_c2_se": s2_se,
            "sigma_c_upper95": float(np.sqrt(max(s2, 0) + 1.96 * s2_se)),
            "note": "elasticities at EZ Rich glucose (c = 0.2) are small, so this bound is weak"}


# ------------------------------------------------------------------ ABC check of (eta2, b)
def abc_check(S, fit, n_draw=4000, seed=1):
    """Simulate y ~ Poisson(q kappa m), m ~ NB(mean M, Fano 1+b) x Gamma(eta2) for the observed genes (exp phase),
    on a grid of (eta2, b); accept the 5% closest (standardised Fano and zero-fraction distances)."""
    rng = np.random.default_rng(seed)
    genes = [(g, d["mean"], d["fano"]) for g, d in fit["genes"].items() if d["mean"] >= 0.01]
    tot = capture_sampler(S, "exp", n_draw, rng); kap = tot / tot.mean()
    obs = {g: (S["operons"][g]["exp"][0], S["operons"][g]["exp"][1], S["operons"][g]["exp"][2][0] / S["cells"]["n_exp"]) for g, _, _ in genes}
    grid = [(e2, b) for e2 in np.linspace(0, 0.4, 17) for b in (0, 0.5, 1, 2, 4, 8)]
    D = []
    for e2, b in grid:
        d = 0.0
        for g, mu, _ in genes:
            M = mu / Q0                                    # molecules per cell at q0
            w = rng.gamma(1 / e2, e2, n_draw) if e2 > 0 else np.ones(n_draw)
            lam = M * w
            m = rng.negative_binomial(lam / b, 1 / (1 + b)) if b > 0 else rng.poisson(lam)
            y = rng.poisson(Q0 * kap * m)
            Fs, z = y.var() / max(y.mean(), 1e-9), (y == 0).mean()
            Fo, zo = obs[g][1] / obs[g][0], obs[g][2]
            nn = S["cells"]["n_exp"]
            sF = np.sqrt((1 / mu + 2 * Fo ** 2) / nn) * np.sqrt(2)
            sz = np.sqrt(zo * (1 - zo) / nn) * np.sqrt(2)
            d += ((Fs - Fo) / sF) ** 2 + ((z - zo) / sz) ** 2
        D.append(d)
    D = np.array(D); keep = D <= np.quantile(D, 0.05)
    acc = np.array(grid)[keep]
    return {"grid_size": len(grid), "accepted": acc.tolist(), "eta2_accepted_range": [float(acc[:, 0].min()), float(acc[:, 0].max())],
            "b_accepted_range": [float(acc[:, 1].min()), float(acc[:, 1].max())],
            "best": {"eta2": float(np.array(grid)[D.argmin()][0]), "b": float(np.array(grid)[D.argmin()][1]), "distance": float(D.min())},
            "n_summary_stats": 2 * len(genes)}


# ------------------------------------------------------------------ 5: lac protein bursts (literature)
def lac_bursts(p):
    """Keep Y_max (and so every A3 deterministic result) fixed; set protein burst k_tl/d_m to Yu 2006 (4.2)."""
    b_old = p["k_tl"] / p["d_m"]
    k_tl = YU2006["burst_size"] * p["d_m"]
    k_tx = p["k_tx"] * p["k_tl"] / k_tl
    q = dict(p, k_tl=k_tl, k_tx=k_tx)
    basal_bursts_per_90min = k_tx / rn.rho_eff(0.0, q) * 90.0      # c = 1, LacI on
    return {"burst_size_before": b_old, "burst_size_after": YU2006["burst_size"], "k_tl": k_tl, "k_tx": k_tx,
            "Y_max_before": rn.y_max(p), "Y_max_after": rn.y_max(q), "basal_LacY": rn.y_max(q) / q["rho"],
            "basal_bursts_per_90min": basal_bursts_per_90min,
            "yu2006": YU2006}


def switch_threshold_lacy(p, T, G=0.0):
    """LacY at the unstable steady state (the hump a cell must cross) at external TMG T: Choi 2008 measured ~375."""
    rho = rn.rho_eff(G, p); n = p["n"]
    a = rn.drive(T, G, p)
    u = np.logspace(-3, 3, 200000)
    f = u * (rho + u ** n) / (1 + u ** n) - a
    roots = u[np.where(np.diff(np.sign(f)))[0]]
    if len(roots) < 3:
        return None
    u_mid = roots[1]
    g = (1 + u_mid ** n) / (rho + u_mid ** n)
    return float(rn.y_max(p) * rn.c_crp(G, p) * g)


# ------------------------------------------------------------------ main
def main(out=ROOT / "outputs"):
    out = Path(out); (out / "a5").mkdir(parents=True, exist_ok=True)
    S = load_summary()
    spec = yaml.safe_load((out / "a3" / "circuit_model.yaml").read_text())
    p = dict(rn.DEFAULT_PARAMS); p.update(spec["parameters"])
    outputs = spec["crp_outputs"]["targets"]
    rng = np.random.default_rng(0)
    R = {"data": str(SUMMARY.relative_to(ROOT)), "conditions": {"G_exp_uM": G_EXP, "G_stat_uM": G_STAT,
         "c_exp": rn.c_crp(G_EXP, p), "c_stat": rn.c_crp(G_STAT, p)}}

    # 1. capture
    cap = {}
    for ph in ("exp", "stat"):
        t = capture_sampler(S, ph, 400000, rng)
        cap[ph] = {"mean_mrna_umi": float(t.mean()), "cv2": float(t.var() / t.mean() ** 2),
                   "relative_to_exp": float(t.mean() / S["mrna_umi"]["exp"]["mean"])}
    R["capture"] = cap
    total_mrna = S["mrna_umi"]["exp"]["mean"] / Q0
    R["molecule_scale"] = {"q0": Q0, "q0_range": Q0_RANGE, "mRNA_per_exp_cell": total_mrna,
                           "mRNA_per_exp_cell_range": [S["mrna_umi"]["exp"]["mean"] / Q0_RANGE[1], S["mrna_umi"]["exp"]["mean"] / Q0_RANGE[0]],
                           "tpm_to_molecules": total_mrna / 1e6, "tpm_to_molecules_A3": spec["crp_outputs"]["tpm_to_molecules"]}
    # 2-3
    R["fano_fit"] = {ph: fit_fano_line(S, ph, cap[ph]["cv2"]) for ph in ("exp", "stat")}
    R["abc"] = abc_check(S, R["fano_fit"]["exp"])
    # 4
    R["crp_noise"] = crp_noise_bound(S, p, outputs)
    # 5
    R["lac_bursts"] = lac_bursts(p)
    q = dict(p, k_tl=R["lac_bursts"]["k_tl"], k_tx=R["lac_bursts"]["k_tx"])
    R["switch_threshold_LacY"] = {f"T{T}": switch_threshold_lacy(q, T) for T in (5, 11, 20, 29)}
    R["switch_threshold_LacY"]["choi2008"] = CHOI2008["switch_threshold_LacY"]
    R["lac_detection"] = {"lacZYA_umi_per_cell": {ph: S["operons"]["lacZYA"][ph][0] for ph in ("exp", "stat")},
                          "fraction_cells_detected": {ph: 1 - S["operons"]["lacZYA"][ph][2][0] / S["cells"]["n_" + ph] for ph in ("exp", "stat")},
                          "model_basal_mRNA": {ph: q["k_tx"] * rn.c_crp(G, q) / q["d_m"] / q["rho"] for ph, G in (("exp", G_EXP), ("stat", G_STAT))},
                          "note": "lac is OFF in PETRI-seq (no inducer): it constrains basal mRNA level, not switching noise"}

    # chosen values
    fe = R["fano_fit"]["exp"]
    eta2 = float(max(fe["eta2"], 0.0))
    chosen = {"output_extrinsic_cv2": eta2, "output_extrinsic_cv2_se": float(fe["eta2_se"]),
              "output_burst_b": 0.0, "output_burst_b_upper95": float(max(fe["burst_b_by_q0"][str(Q0)]["b_upper95"], 0.0)),
              "crp_extrinsic_sigma": 0.0, "crp_extrinsic_sigma_upper95": R["crp_noise"]["sigma_c_upper95"],
              "lac_protein_burst": YU2006["burst_size"], "tpm_to_molecules": R["molecule_scale"]["tpm_to_molecules"],
              "capture_per_molecule_q0": Q0}
    R["chosen"] = chosen

    # calibrated YAML: A3 v3 + new parameters + noise + observation model
    cal = copy.deepcopy(spec)
    cal["name"] = "crp_lac_circuit_A5_calibrated"
    cal["version"] = "A5 calibrated (24 Nov 2026), on A3 v3"
    cal["description"] = ("A3 v3 circuit with noise parameters calibrated in A5 (PETRI-seq GSE141018 exp 2.01 + Yu 2006). "
                          "Deterministic behaviour unchanged (Y_max, rho, thresholds identical); noise changed. "
                          "Use this file, not circuit_model.yaml, from A5 on (A6, Student B Module 4).")
    cal["parameters"]["k_tl"] = float(R["lac_bursts"]["k_tl"])
    cal["parameters"]["k_tx"] = float(R["lac_bursts"]["k_tx"])
    cal["crp_outputs"]["tpm_to_molecules"] = float(chosen["tpm_to_molecules"])
    cal["noise"] = {
        "parameters": {k: chosen[k] for k in ("output_extrinsic_cv2", "output_burst_b", "crp_extrinsic_sigma", "lac_protein_burst")},
        "uncertainty": {"output_extrinsic_cv2_se": chosen["output_extrinsic_cv2_se"], "output_burst_b_upper95": chosen["output_burst_b_upper95"],
                        "crp_extrinsic_sigma_upper95": chosen["crp_extrinsic_sigma_upper95"], "lac_protein_burst_se": YU2006["burst_size_se"]},
        "sources": {"output_extrinsic_cv2": "PETRI-seq exp 2.01: slope of Fano-1 vs mean after capture (method of moments; ABC check)",
                    "output_burst_b": "PETRI-seq: intercept not distinguishable from Poisson at ~5% capture (upper bound given)",
                    "crp_extrinsic_sigma": "PETRI-seq: CRP targets uncorrelated after capture scaling (upper bound given)",
                    "lac_protein_burst": "Yu et al. 2006 Science 311:1600 (repressed lac promoter, 4.2 +/- 0.5 per burst)",
                    "k_tl, k_tx": "k_tl = 4.2 d_m; k_tx raised so Y_max = k_tx k_tl / (d_m d_y) is unchanged",
                    "tpm_to_molecules": f"PETRI-seq mean {S['mrna_umi']['exp']['mean']:.0f} mRNA UMIs per exp cell / capture q0 = {Q0}"},
        "changed_from_A3": {"k_tl": [p["k_tl"], R["lac_bursts"]["k_tl"]], "k_tx": [p["k_tx"], R["lac_bursts"]["k_tx"]],
                            "tpm_to_molecules": [spec["crp_outputs"]["tpm_to_molecules"], chosen["tpm_to_molecules"]]},
    }
    cal["observation_model"] = {
        "description": "To compare simulated cells with PETRI-seq: UMI_ig ~ Poisson(q_g * kappa_i * molecules_ig); "
                       "kappa_i = cell's total mRNA UMIs / exp-phase mean, drawn from the phase's empirical distribution.",
        "capture": {ph: {**cap[ph], "hist_log2_halfbins": S["mrna_umi"][ph]["hist_log2_halfbins"]} for ph in cap},
        "q0": Q0, "q_g": "per operon, fitted on exponential-phase means in a5_validate.py (outputs/a5/validation.json)",
    }
    cal = json.loads(json.dumps(cal, default=float))
    (out / "a5" / "calibrated_circuit_model.yaml").write_text(yaml.safe_dump(cal, sort_keys=False, allow_unicode=True, width=120))
    (out / "a5" / "noise_fit.json").write_text(json.dumps(R, indent=2, default=float))
    return R


if __name__ == "__main__":
    R = main(*sys.argv[1:])
    print(json.dumps({k: R[k] for k in ("capture", "molecule_scale", "chosen", "switch_threshold_LacY", "lac_bursts")}, indent=1, default=float))
    fe = R["fano_fit"]["exp"]; print("eta2", fe["eta2"], "+/-", fe["eta2_se"], "A", fe["intercept_A"], "ABC", R["abc"]["best"], R["abc"]["eta2_accepted_range"], R["abc"]["b_accepted_range"])
    print("crp sigma_c2", R["crp_noise"]["sigma_c2_hat"], "+/-", R["crp_noise"]["sigma_c2_se"], "upper", R["crp_noise"]["sigma_c_upper95"])
