"""A3 regression tests: python -m pytest tests/"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import reaction_network as rn  # noqa: E402

P = json.loads((ROOT / "outputs" / "a3" / "fitted_params.json").read_text())["params"]


def test_rho_equals_max_over_basal():
    up = rn.simulate(1e-3, 0, P, induced=False, t_end=200 * 60)
    full = rn.y_max(P)
    assert full / up == pytest.approx(P["rho"], rel=0.02)


def test_cusp_at_nine():
    assert rn.switching_drives(8.9) is None
    assert rn.switching_drives(9.2) is not None


def test_ode_converges_to_analytic_thresholds():
    off, on = rn.thresholds_tmg(0.0, P)
    t = 300 * 60
    assert rn.simulate(on * 1.05, 0, P, induced=False, t_end=t) > 100      # switches on just above
    assert rn.simulate(on * 0.95, 0, P, induced=False, t_end=t) < 100      # stays off just below
    assert rn.simulate(off * 1.05, 0, P, induced=True, t_end=t) > 100      # stays on just above
    assert rn.simulate(off * 0.95, 0, P, induced=True, t_end=t) < 100      # switches off just below


def test_bistable_band_contains_two_states():
    T = 12.0
    lo = rn.simulate(T, 0, P, induced=False)
    hi = rn.simulate(T, 0, P, induced=True)
    assert hi / lo > 50


def test_glucose_raises_thresholds():
    off0, on0 = rn.thresholds_tmg(0, P)
    off1, on1 = rn.thresholds_tmg(1000, P)
    assert on1 > 10 * on0 and off1 > 10 * off0


def test_rho_falls_with_crp():
    """Kuhlman 2007: CRP strengthens LacI looping, so glucose (less CRP) and crp knockout lower rho."""
    r0, r1, rko = rn.rho_eff(0, P), rn.rho_eff(1000, P), rn.rho_eff(0, P, crp_knockout=True)
    assert P["rho_crp_coupling"] == "knockout"
    assert r0 == r1 == pytest.approx(P["rho"])      # glucose leaves rho alone (keeps the Fig 2c fit)
    assert r0 > 5 * rko and rko > 9                 # crp knockout: weaker loop, still a switch
    g = dict(P, rho_crp_coupling="glucose")
    assert rn.rho_eff(0, g) > rn.rho_eff(1000, g) > rn.rho_eff(0, g, crp_knockout=True)
    # full CRP range spans Kuhlman's ~10x change in rho - 1 (Omega = 12.4, L0 = 3.9)
    full = 1 + P["loop_L0"] * P["loop_omega"]
    assert (full / (1 + P["loop_L0"])) == pytest.approx(10.0, rel=0.05)


def test_fixed_rho_switch_restores_old_model():
    q = dict(P, rho_crp_coupling=False)
    assert rn.rho_eff(1000, q) == q["rho"]
    assert rn.k_rep_max(q, 1000) == pytest.approx((q["rho"] - 1) * q["k_unrep"])


def test_gillespy2_model_builds():
    gp = pytest.importorskip("gillespy2")
    m = rn.to_gillespy2(P, 15.0, 0.0)
    assert len(m.listOfReactions) == 8


# ---------------- A3 v3 additions
def test_crp_knockout_stays_off():
    """Kuhlman 2007 / Traxler 2006: without CRP lac does not induce at any normal TMG dose."""
    for T in [10, 100, 1000]:
        assert rn.simulate(T, 0, P, induced=False, crp_knockout=True) < 5


def test_crp_free_basal_from_kuhlman():
    assert rn.c_eff(0, P, crp_knockout=True) == pytest.approx(1 / 240)
    assert rn.c_eff(0, P) == pytest.approx(1.0)


def test_sbml_loads_in_gillespy2():
    import gillespy2
    sbml = ROOT / "outputs" / "a3" / "circuit_model.xml"
    m, errs = gillespy2.import_SBML(str(sbml))
    assert not errs and {"P_free", "P_rep", "M", "Y", "X"} <= set(m.listOfSpecies)
    assert len(m.listOfReactions) == 8 + 2 * 21
