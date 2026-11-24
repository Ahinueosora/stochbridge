"""A2: the mediation rule on synthetic data (error rates) and, if PRECISE-1K is present, the published 21 edges."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))
from mediation_test import ols  # noqa: E402

A2_EDGES = {"hflC", "tdcB", "uhpT", "cas3", "rbsC", "yahN", "fecE", "mhpF", "flhC", "flhD", "araJ", "rbsA", "nadC",
            "galK", "rbsB", "gyrA", "rbsR", "lamB", "malK", "preA", "idnO"}


def _rule(y, D, C):
    one = np.ones((len(y), 1))
    b, p = ols(y, np.hstack([one, D[:, None]]))
    bc, pc = ols(y, np.hstack([one, D[:, None], C[:, None]]))
    return p[1] < 0.05 and pc[1] > 0.05 and pc[2] < 0.05 and (b[1] - bc[1]) / b[1] >= 0.5


def test_rule_accepts_mediated_and_rejects_direct_effects():
    rng = np.random.default_rng(0)
    D = np.r_[np.zeros(3), np.ones(17)]
    hits_med = hits_dir = 0
    for _ in range(400):
        C = -2 * D + rng.normal(0, 0.6, 20)
        hits_med += _rule(1.0 * C + rng.normal(0, 0.3, 20), D, C)        # effect only through C
        hits_dir += _rule(-2 * D + rng.normal(0, 0.3, 20), D, C)         # effect bypasses C
    assert hits_med / 400 > 0.8 and hits_dir / 400 < 0.05


@pytest.mark.skipif(not (ROOT / "outputs" / "a2" / "mediation_summary.json").exists(), reason="run src/mediation_test.py")
def test_real_data_gives_the_21_edges():
    s = json.loads((ROOT / "outputs" / "a2" / "mediation_summary.json").read_text())
    assert set(s["causal_genes"]) == A2_EDGES
    assert s["perturbation_effect"] == 188 and s["causal_with_mRNA_mediator"] == 0
    assert s["r_leave_none_out_vs_published_Crp1"] > 0.999
