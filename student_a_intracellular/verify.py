"""
Verify Student A's results (StochBridge Track A, tasks A1-A6).

    python get_data.py      # once: downloads PRECISE-1K (pinned version, checksums checked)
    python verify.py        # about 5 minutes

Three kinds of check, each printed as PASS / FAIL:
  1. RERUN   A1 and A2 are recomputed from the raw data and compared with the files in outputs/
             (A3 too with --a3, adds about 3 minutes).
  2. TESTS   the pytest suite (38 tests: circuit, fixed-seed simulator regression, noise layer, perturbation API,
             and the Lamoureux Fig. 6 reproduction).
  3. RESULTS each task's headline numbers, read from outputs/, compared with the numbers quoted in the report.

To regenerate every output from scratch instead (about 25 min): python run_all.py,
then python verify.py.

Options:  --a3       also recompute A3 (lac circuit fit + final model)
          --no-tests skip pytest
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
PY = sys.executable
RESULTS = []


def check(task, what, ok, detail=""):
    RESULTS.append((task, what, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {task:3} {what}" + (f"  [{detail}]" if detail else ""), flush=True)


def load(rel):
    p = OUT / rel
    return json.loads(p.read_text()) if p.suffix == ".json" else p


def close(a, b, tol):
    return a is not None and b is not None and abs(a - b) <= tol


# ---------------------------------------------------------------- 0. data
def check_data():
    print("\n0. DATA")
    import get_data
    for rel, want in get_data.FILES.items():
        p = get_data.DEST / rel
        if not p.exists():
            check("-", f"data file {rel}", False, "missing: run `python get_data.py`")
            sys.exit(1)
    bad = [rel for rel, want in get_data.FILES.items() if get_data.sha256(get_data.DEST / rel) != want]
    check("-", f"{len(get_data.FILES)} PRECISE-1K files match pinned SHA-256 checksums", not bad, ", ".join(bad))
    petri = list((ROOT / "data/scrnaseq/GSE141018").glob("GSM4489500_*.txt.gz"))
    check("-", "PETRI-seq file present (GEO GSM4489500)", petri)


# ---------------------------------------------------------------- 1. rerun
def rerun(with_a3):
    import pandas as pd
    print("\n1. RERUN FROM RAW DATA (into a temporary copy; outputs/ is not touched)")
    tmp = Path(tempfile.mkdtemp(prefix="verify_"))
    try:
        for d in ["src", "notebooks", "config"]:
            shutil.copytree(ROOT / d, tmp / d)
        (tmp / "data").mkdir()
        for d in ["precise1k", "scrnaseq", "curated"]:
            (tmp / "data" / d).symlink_to(ROOT / "data" / d, target_is_directory=True)
        for a in ["a1", "a2", "a3"]:
            (tmp / "outputs" / a).mkdir(parents=True)
        import os
        env = {**os.environ, "MPLBACKEND": "Agg", "PATH": os.path.dirname(PY) + os.pathsep + os.environ.get("PATH", "")}

        def run(cmd, cwd=tmp):
            r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
            if r.returncode:
                print(r.stdout[-2000:], r.stderr[-2000:])
            return r.returncode == 0

        ok = run([PY, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--output",
                  str(tmp / "a1_check.ipynb"), "notebooks/01_crp_reproduction.ipynb"])
        new = pd.read_csv(tmp / "outputs/a1/crp_activity_by_condition.csv", index_col=0) if ok else None
        old = pd.read_csv(OUT / "a1/crp_activity_by_condition.csv", index_col=0)
        diff = float((new[["Crp-1", "Crp-2"]] - old.loc[new.index, ["Crp-1", "Crp-2"]]).abs().max().max()) if ok else None
        check("A1", "Crp-1/Crp-2 activities of all 1,035 samples recomputed = committed", ok and diff < 1e-9,
              f"max difference {diff:.1e}" if ok else "notebook failed")
        if ok:
            counts = new.crp_category.value_counts().to_dict()
            check("A1", "CRP mutant groups: 9 deletion, 3 ΔAr1+Ar2, 3 ΔAr2, 2 ΔAr1",
                  counts.get("CRP full deletion") == 9 and counts.get("CRP ΔAr1+Ar2") == 3
                  and counts.get("CRP ΔAr2 (Class II)") == 3 and counts.get("CRP ΔAr1 (Class I)") == 2, str(counts))

        ok = run([PY, "src/mediation_test.py"])
        new = json.loads((tmp / "outputs/a2/mediation_summary.json").read_text()) if ok else {}
        check("A2", "causal-mediation test recomputed: summary identical to committed",
              ok and new == load("a2/mediation_summary.json"),
              f"{new.get('causal')} causal edges of {new.get('targets_tested')} genes tested" if ok else "script failed")
        if ok:
            a = pd.read_csv(tmp / "outputs/a2/mediation_results.csv")
            b = pd.read_csv(OUT / "a2/mediation_results.csv")
            check("A2", "per-gene statistics (mediation_results.csv) identical", a.equals(b))

        if with_a3:
            ok = run([PY, "src/run_a3.py", str(ROOT / "data/precise1k/precise1k/data"), "outputs/a3"]) \
                and run([PY, "../../src/build_network.py", "../../data/curated/lac_network_gene_check.xlsx"],
                        cwd=tmp / "outputs/a3") \
                and run([PY, "src/build_a3_v3.py", "outputs/a3/network_v3_edges.csv", "outputs/a3"])
            if not ok:
                check("A3", "A3 rerun", False, "script failed")
            else:
                import yaml
                n = yaml.safe_load((tmp / "outputs/a3/circuit_model.yaml").read_text())["parameters"]
                o = yaml.safe_load((OUT / "a3/circuit_model.yaml").read_text())["parameters"]
                num = [k for k in o if isinstance(o[k], (int, float)) and not isinstance(o[k], bool)]
                worst = max(abs(n[k] - o[k]) / max(abs(o[k]), 1e-12) for k in num)
                check("A3", f"lac circuit parameters recomputed ({len(num)} values) = committed",
                      worst < 1e-5, f"largest relative difference {worst:.1e} (optimiser tolerance)")
                a = json.loads((tmp / "outputs/a3/validated_edges.json").read_text())
                check("A3", "circuit edge list with evidence identical", a == load("a3/validated_edges.json"))
                f = pd.read_csv(tmp / "outputs/a3/crp_outputs_fit.csv")
                g = pd.read_csv(OUT / "a3/crp_outputs_fit.csv")
                cols = f.select_dtypes("number").columns
                worst = float(((f[cols] - g[cols]).abs() / g[cols].abs().clip(lower=1e-12)).max().max())
                check("A3", "21 CRP-output gene fits recomputed = committed within 1%", worst < 0.01,
                      f"largest relative difference {worst:.1e}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- 2. tests
def tests():
    print("\n2. TESTS (pytest, about 2 minutes)")
    r = subprocess.run([PY, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests"], cwd=ROOT,
                       capture_output=True, text=True)
    last = [l for l in r.stdout.strip().splitlines() if l.strip()][-1]
    check("-", "pytest suite", r.returncode == 0, last)


# ---------------------------------------------------------------- 3. results
def results():
    import pandas as pd
    import yaml
    print("\n3. RESULTS vs the numbers in the report")

    s = load("a2/mediation_summary.json")
    check("A2", "21 causal CRP edges out of 543 genes tested; 188 respond to CRP loss",
          (s["causal"], s["targets_tested"], s["perturbation_effect"]) == (21, 543, 188))
    check("A2", "leave-none-out CRP activity reproduces published Crp-1 (r > 0.999)",
          s["r_leave_none_out_vs_published_Crp1"] > 0.999, f"r = {s['r_leave_none_out_vs_published_Crp1']:.4f}")

    v = yaml.safe_load((OUT / "a3/circuit_model.yaml").read_text())["validation"]
    off = v["ozbudak_fig2c"]["switch_off_rmse_log10_heldout"]
    check("A3", "held-out Ozbudak switch-OFF thresholds predicted within ~23% (RMSE of log10)",
          10 ** off - 1 < 0.25, f"{100 * (10 ** off - 1):.0f}%")
    c = load("a3/a3_checks_v3.json")
    check("A3", "bistability cusp at rho ≈ 9", 8 <= c["cusp_rho"] <= 10, f"rho = {c['cusp_rho']:.2f}")
    check("A3", "rho = 170 and 50 strains hysteretic, rho = 5 strain not",
          c["strain_rho170"]["hysteretic"] and c["strain_rho50"]["hysteretic"] and not c["strain_rho5"]["hysteretic"])
    check("A3", "crp knockout stays OFF (max LacY < 1% of wild type)",
          c["crp_knockout"]["max_LacY_any_TMG_up_to_2000uM"] < 0.01 * c["crp_knockout"]["wt_max_LacY"])

    a4 = load("a4/vilar_results.json")
    h = a4["ssa_hysteresis"]
    check("A4", "SSA switch-on 37 uM and switch-off 3.4 uM (Ozbudak: 30 and 3)",
          close(h["switch_on_T50_uM"], 37.3, 0.5) and close(h["switch_off_T50_uM"], 3.41, 0.05),
          f"{h['switch_on_T50_uM']:.1f} / {h['switch_off_T50_uM']:.2f} uM")
    b = a4["bimodality_wt"]["1200"]
    check("A4", "normal switch (rho 170) bimodal at 20 h (dip test p < 0.01)", b["dip_p"] < 0.01,
          f"p = {b['dip_p']:.3g}, {100 * b['fraction_ON']:.0f}% ON")
    r5 = a4["rho5_control"]["bimodality_rho5"]
    check("A4", "switch-less strain (rho 5) stays unimodal at every time (dip p > 0.05)",
          all(x["dip_p"] > 0.05 for x in r5.values()))
    check("A4", "maintenance dose 11 uM: 17.4% of OFF cells switch, 0% of ON cells switch back",
          close(a4["vilar"]["vilar3b_maint_from_OFF"]["fraction_switched_30h"], 0.174, 0.001)
          and a4["vilar"]["vilar3b_maint_from_ON"]["fraction_switched_30h"] == 0.0)
    fe = load("a4/fast_vs_exact.json")
    z = [abs(fe[k]["switched_fraction_diff"]) / fe[k]["switched_fraction_diff_se"]
         for k in fe if isinstance(fe[k], dict) and fe[k].get("switched_fraction_diff_se")]
    check("A4", "fast mode matches exact mode (switched fraction within 2 SE)", all(x < 2 for x in z),
          ", ".join(f"{x:.1f} SE" for x in z))

    a5 = load("a5/validation.json")
    check("A5", "real PETRI-seq: 17-operon PC1 separates growth phases (AUC 0.81)",
          close(a5["real"]["pca17_auc"], 0.806, 0.005), f"AUC = {a5['real']['pca17_auc']}")
    m = a5["memory_check"]
    check("A5", "calibration cuts false lac switching from 19% to under 1%",
          m["uncalibrated_A4"]["uninduced"]["fraction_flipped_30h"] >= 0.15
          and m["calibrated"]["uninduced"]["fraction_flipped_30h"] < 0.01,
          f"{100 * m['uncalibrated_A4']['uninduced']['fraction_flipped_30h']:.0f}% -> "
          f"{100 * m['calibrated']['uninduced']['fraction_flipped_30h']:.1f}%")
    ag = a5["calibrated"]["agreement"]["exp"]["log10_mean_r"]
    check("A5", "calibrated means match growing cells (r of log10 means > 0.99)", ag > 0.99, f"r = {ag:.4f}")

    aat = load("a6/aat_check.json")
    check("A6", "Lamoureux Fig. 6 reproduced: all 11 paper DiMAs found, max error < 0.02",
          len(aat["overlap"]) == 11 and not aat["only_paper"] and aat["max_abs_diff_vs_paper"] < 0.02,
          f"{len(aat['overlap'])}/11, max error {aat['max_abs_diff_vs_paper']:.3f}")
    S = pd.read_csv(OUT / "a6/perturbation_summary.csv")
    rec, sim, real = S.panel_recovered.sum(), S.panel_sim_dimas.sum(), S.panel_real_dimas.sum()
    check("A6", "circuit-panel DiMA: 16 of 22 simulated calls correct, 16 of 48 real changes found",
          (rec, sim, real) == (16, 22, 48), f"{rec}/{sim} correct, {rec}/{real} found")
    check("A6", "stricter Fig. 6 thresholds: 8 of 9 correct; multi-gene iModulons only: 10 of 12",
          (S.panel_fig6_recovered.sum(), S.panel_fig6_sim_dimas.sum()) == (8, 9)
          and (S.panel_multigene_recovered.sum(), S.panel_multigene_sim_dimas.sum()) == (10, 12))
    H = pd.read_csv(OUT / "a6/heldout/perturbation_summary.csv")
    P = S[S.perturbation.isin(H.perturbation)]
    check("A6", "held-out refit (CRP curves refitted without the 32 compared samples): 9 of 11 correct vs 10 of 12",
          (H.panel_recovered.sum(), H.panel_sim_dimas.sum(), P.panel_recovered.sum(), P.panel_sim_dimas.sum())
          == (9, 11, 10, 12),
          f"held out {H.panel_recovered.sum()}/{H.panel_sim_dimas.sum()}, in-sample {P.panel_recovered.sum()}/{P.panel_sim_dimas.sum()}")


if __name__ == "__main__":
    sys.path.insert(0, str(ROOT))
    check_data()
    rerun(with_a3="--a3" in sys.argv)
    if "--no-tests" not in sys.argv:
        tests()
    results()
    n_fail = sum(not ok for *_, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail} of {len(RESULTS)} checks passed." + ("" if not n_fail else "  See FAIL lines above."))
    sys.exit(1 if n_fail else 0)
