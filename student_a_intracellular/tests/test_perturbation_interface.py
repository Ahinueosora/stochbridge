"""A6 / Module 2: the Perturbation API contract Student B's screen_runner relies on."""
import copy
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))
from perturbation import (CarbonSourceShift, GeneKnockout, InducerLevelChange, Perturbation,  # noqa: E402
                          load_condition_map, load_perturbations)
from ssa_simulator import StochasticCircuit  # noqa: E402

MODEL = yaml.safe_load((ROOT / "outputs" / "calibrated_circuit_model.yaml").read_text())


def test_apply_returns_new_model_and_leaves_input_alone():
    before = copy.deepcopy(MODEL)
    m = GeneKnockout("crp").apply(MODEL)
    assert MODEL == before
    assert m["genetics"]["crp_knockout"] and m["genetics"]["deleted_genes"] == ["crp"]
    assert m["perturbation"]["kind"] == "GeneKnockout" and m["name"].endswith("__dcrp")
    assert Perturbation().apply(MODEL)["parameters"] == MODEL["parameters"]


def test_knockout_effects_and_refusals():
    assert GeneKnockout("lacI").apply(MODEL)["parameters"]["rho"] == 1.0
    assert GeneKnockout("lacY").apply(MODEL)["parameters"]["k_tl"] == 0.0
    t = GeneKnockout("malK").apply(MODEL)["crp_outputs"]["targets"]["malK"]
    assert t["basal_TPM"] == 0 and t["amplitude_TPM"] == 0
    with pytest.raises(NotImplementedError):
        GeneKnockout("lacZ")
    with pytest.raises(NotImplementedError):
        CarbonSourceShift(carbon_source="fructose")
    assert abs(CarbonSourceShift(carbon_source="glucose", g_per_L=2).glucose_uM - 11101.9) < 1
    assert CarbonSourceShift(carbon_source="glycerol").glucose_uM == 0


def test_applied_model_runs_same_as_direct_call():
    p = GeneKnockout("crp", n_cells=40, t_end_min=300, dt_min=100, seed=4)
    sc = StochasticCircuit(MODEL)
    a = StochasticCircuit(p.apply(MODEL)).run().data
    b = sc.run(p).data
    assert a.equals(b)
    assert not sc.run(p, seed=5).data.equals(b)          # replicate populations differ by seed only


def test_iptg_induces_without_lacY_feedback():
    sc = StochasticCircuit(MODEL)
    G = CarbonSourceShift(carbon_source="glucose", g_per_L=2).glucose_uM
    on = sc.run(Perturbation(IPTG_uM=50, glucose_uM=G, n_cells=100, seed=2))
    y_on = on.meta["Y_on"]
    assert on.at(1200).Y_count.mean() > 0.8 * y_on
    ko = sc.run(Perturbation(IPTG_uM=50, glucose_uM=G, n_cells=50, seed=2, param_overrides={"k_tl": 0.0}))
    assert ko.at(1200).Y_count.max() == 0


def test_condition_map_loads_every_perturbation_with_a_control():
    cmap = load_condition_map()
    trip = load_perturbations(with_controls=True)
    assert len(trip) == len(cmap["perturbations"]) >= 8
    names = [p.name for p, _, _ in trip]
    assert len(set(names)) == len(names)
    for p, c, m in trip:
        assert c.name == p.name + "__control"
        if m:
            assert {"dataset", "control", "perturbed", "held_out"} <= set(m)


@pytest.mark.skipif(not (ROOT / "data" / "precise1k" / "p1k_metadata_qc.csv").exists(), reason="PRECISE-1K not staged")
def test_every_matched_condition_exists_in_the_data():
    p1k = pd.read_csv(ROOT / "data" / "precise1k" / "p1k_metadata_qc.csv", index_col=0)
    k12 = pd.read_csv(ROOT / "data" / "precise1k" / "k12_metadata_qc.csv", index_col=0, low_memory=False)
    for _, _, m in load_perturbations(with_controls=True):
        if not m:
            continue
        md = p1k if m["dataset"] == "precise1k" else k12
        for key in (m["control"], m["perturbed"]):
            proj, cond = key.split(":")
            assert ((md.project == proj) & (md.condition == cond)).sum() >= 2, key


def test_bh_matches_definition():
    from dima import bh
    p = pd.Series([0.01, 0.04, 0.03, 0.2])
    q = bh(p)
    assert np.allclose(q.values, [0.04, 0.0533333, 0.0533333, 0.2])


@pytest.mark.skipif(not (ROOT / "outputs" / "a6" / "aat_check.json").exists(), reason="run src/a6_reproduce_aat.py")
def test_fig6_aat_reproduced():
    c = json.loads((ROOT / "outputs" / "a6" / "aat_check.json").read_text())
    assert len(c["overlap"]) == 11 and not c["only_paper"] and c["max_abs_diff_vs_paper"] < 0.05


def test_subclasses_round_trip_through_to_dict():
    for p in (InducerLevelChange(40, n_cells=5, seed=3), InducerLevelChange(IPTG_uM=50), CarbonSourceShift(carbon_source="glycerol"),
              GeneKnockout("cyaA")):
        d = p.to_dict()
        assert type(p)(**d).to_dict() == d
        assert Perturbation.from_dict(d).to_dict() == d
