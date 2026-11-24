"""A4 / Module 1 regression tests: python -m pytest -q tests/test_ssa_reproducibility.py

Fixed seeds must give bit-identical trajectories; a stored golden summary catches silent changes to the model,
the engine or GillesPy2. Regenerate the golden file only on purpose:  python tests/test_ssa_reproducibility.py --update-golden
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from io_schema import CORE_COLUMNS, Trajectory, validate  # noqa: E402
from perturbation import GeneKnockout, InducerLevelChange, OperatorDecoy, Perturbation  # noqa: E402
from ssa_simulator import StochasticCircuit  # noqa: E402

YAML = ROOT / "outputs" / "circuit_model.yaml"
GOLDEN = Path(__file__).parent / "golden_ssa_fast.json"
GOLDEN_CASES = [InducerLevelChange(40, n_cells=50, t_end_min=600, seed=2026, name="g_T40"),
                InducerLevelChange(11, initial_state="induced", n_cells=30, t_end_min=600, seed=7, name="g_T11_ON"),
                GeneKnockout("crp", TMG_uM=100, n_cells=30, t_end_min=600, seed=3, name="g_dcrp")]


@pytest.fixture(scope="module")
def sc():
    return StochasticCircuit(YAML, mode="fast", include_outputs=True)


def _fingerprint(tr):
    d = tr.data
    return {"rows": int(len(d)), "Y_sum": int(d.Y_count.sum()), "M_sum": int(d.M_count.sum()),
            "P_free_sum": int(d.P_free_count.sum()), "final_Y": d[d.time == d.time.max()].Y_count.tolist()[:10]}


def test_same_seed_identical(sc):
    p = InducerLevelChange(40, n_cells=20, t_end_min=600, seed=11)
    a, b = sc.run(p), sc.run(p)
    assert a.data.equals(b.data)


def test_different_seed_differs(sc):
    a = sc.run(InducerLevelChange(40, n_cells=20, t_end_min=600, seed=11))
    b = sc.run(InducerLevelChange(40, n_cells=20, t_end_min=600, seed=12))
    assert not a.data.equals(b.data)


def test_golden_regression(sc):
    if not GOLDEN.exists():
        pytest.skip("no golden file; run with --update-golden")
    gold = json.loads(GOLDEN.read_text())
    for p in GOLDEN_CASES:
        assert _fingerprint(sc.run(p)) == gold[p.name], p.name


def test_schema_and_parquet_roundtrip(sc, tmp_path):
    tr = sc.run(InducerLevelChange(15, n_cells=5, t_end_min=120, seed=1))
    assert list(tr.data.columns[:len(CORE_COLUMNS)]) == CORE_COLUMNS
    assert len(tr.data.columns) == len(CORE_COLUMNS) + 21
    validate(tr.data)
    f = tr.to_parquet(tmp_path / "t.parquet")
    back = Trajectory.from_parquet(f)
    assert back.data.equals(tr.data)
    assert back.meta["perturbation"] == tr.meta["perturbation"]
    assert back.meta["schema_version"].startswith("1.")
    assert set(tr.times) == set(np.arange(0, 121, 5.0))


def test_schema_rejects_bad_frames(sc):
    tr = sc.run(InducerLevelChange(15, n_cells=2, t_end_min=20, seed=1))
    with pytest.raises(ValueError):
        validate(tr.data.rename(columns={"Y_count": "LacY"}))
    with pytest.raises(ValueError):
        validate(tr.data.astype({"Y_count": "float64"}))


# ---------------- behaviour the SSA must show
def test_basal_mean_matches_ode(sc):
    tr = sc.run(InducerLevelChange(1, n_cells=400, t_end_min=1200, seed=5))
    y = tr.at(1200).Y_count.values
    target = tr.meta["Y_off"]                         # deterministic basal LacY
    assert abs(y.mean() - target) < 4 * y.std(ddof=1) / np.sqrt(len(y)) + 0.5


def test_all_cells_on_at_high_tmg(sc):
    tr = sc.run(InducerLevelChange(500, n_cells=200, t_end_min=1800, seed=6))
    assert (tr.at(1800).Y_count.values > tr.meta["on_threshold_Y"]).mean() > 0.95


def test_memory_inside_zone(sc):
    """Ozbudak: cells that start ON stay ON at a dose where cells that start OFF stay mostly OFF."""
    on = sc.run(InducerLevelChange(11, initial_state="induced", n_cells=200, seed=8))
    off = sc.run(InducerLevelChange(11, initial_state="uninduced", n_cells=200, seed=9))
    thr = on.meta["on_threshold_Y"]
    assert (on.at(1200).Y_count.values > thr).mean() > 0.9
    assert (off.at(1200).Y_count.values > thr).mean() < 0.3


def test_bimodal_not_unimodal(sc):
    """Wild type at a switching dose: two peaks (dip test). Graded rho = 5 strain: one peak."""
    from io_schema import bimodality_of_counts
    wt = sc.run(InducerLevelChange(40, n_cells=600, t_end_min=1200, seed=10))
    g = sc.run(OperatorDecoy(5, TMG_uM=3, n_cells=600, t_end_min=1200, seed=10))
    assert bimodality_of_counts(wt.at(1200).Y_count.values)["dip_p"] < 0.01
    assert bimodality_of_counts(g.at(1200).Y_count.values)["dip_p"] > 0.05


def test_crp_knockout_stays_off(sc):
    tr = sc.run(GeneKnockout("crp", TMG_uM=500, n_cells=100, seed=4))
    assert tr.at(1200).Y_count.mean() < 2


def test_perturbation_validation():
    with pytest.raises(ValueError):
        Perturbation(initial_state="sideways")
    with pytest.raises(NotImplementedError):
        GeneKnockout("lacZ")


def test_exact_mode_builds_and_is_reproducible():
    ex = StochasticCircuit(YAML, mode="exact", include_outputs=False)
    p = InducerLevelChange(1, n_cells=3, t_end_min=120, seed=3)
    a, b = ex.run(p), ex.run(p)
    assert a.data.equals(b.data) and a.meta["x_column"] == "simulated"


if __name__ == "__main__" and "--update-golden" in sys.argv:
    s = StochasticCircuit(YAML, mode="fast", include_outputs=True)
    GOLDEN.write_text(json.dumps({p.name: _fingerprint(s.run(p)) for p in GOLDEN_CASES}, indent=1))
    print("wrote", GOLDEN)
