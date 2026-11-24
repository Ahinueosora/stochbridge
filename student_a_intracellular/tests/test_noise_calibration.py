"""A5: calibrated model loads, noise layer is used, sampler is reproducible and has the intended moments."""
from pathlib import Path
import sys

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))
from perturbation import Perturbation  # noqa: E402
from ssa_simulator import StochasticCircuit  # noqa: E402

CAL = ROOT / "outputs" / "calibrated_circuit_model.yaml"


def test_yaml_has_noise_layer_on_top_of_a3():
    y = yaml.safe_load(CAL.read_text())
    a3 = yaml.safe_load((ROOT / "outputs" / "circuit_model.yaml").read_text())
    n = y["noise"]["parameters"]
    assert 0 < n["output_extrinsic_cv2"] < 0.5 and n["lac_protein_burst"] == 4.2
    assert len(y["circuit_edges"] if "circuit_edges" in y else y["edges"]) == len(a3["circuit_edges"] if "circuit_edges" in a3 else a3["edges"])
    assert "petri_ez_rich_medium" in y["observation_model"]


def test_outputs_sampler_reproducible_and_moments():
    sc = StochasticCircuit(CAL, mode="fast", include_outputs=True)
    pt = Perturbation(TMG_uM=0.0, glucose_uM=11100, n_cells=3000, t_end_min=300, dt_min=300, seed=5)
    a, b = sc.run(pt), sc.run(pt)
    col = [c for c in a.data.columns if c.startswith("mRNA_")][0]
    assert np.array_equal(a.data[col].values, b.data[col].values)
    assert a.meta["crp_outputs_method"].startswith("exact")
    x = a.at(300)[col].values.astype(float)
    m = x.mean()
    eta2 = sc.noise["output_extrinsic_cv2"]
    fano = x.var() / m
    assert abs(fano - (1 + eta2 * m)) < 0.25 * (1 + eta2 * m) + 0.1   # negative-binomial-like spread


def test_lac_leak_reduced():
    sc = StochasticCircuit(CAL, mode="fast", include_outputs=False)
    from perturbation import InducerLevelChange
    tr = sc.run(InducerLevelChange(11.0, initial_state="uninduced", n_cells=300, t_end_min=1800, dt_min=30, seed=3))
    W = tr.wide("Y").values
    assert (W > tr.meta["on_threshold_Y"]).any(1).mean() < 0.05


RAW = ROOT / "data" / "scrnaseq" / "GSE141018"


@pytest.mark.skipif(not any(RAW.glob("GSM4489500*.txt.gz")), reason="raw GEO file not present")
def test_python_loader_reproduces_the_summary_used_for_calibration():
    import json
    from load_scrnaseq import load_petri, summarise
    raw = next(RAW.glob("GSM4489500*.txt.gz"))
    cells, X, operons = load_petri(raw)
    b = summarise(cells, X, operons)
    a = json.loads((RAW / "petri_2_01_summary.json").read_text())
    assert (b["cells"]["n_exp"], b["cells"]["n_stat"]) == (a["cells"]["n_exp"], a["cells"]["n_stat"]) == (3430, 4499)
    for g in b["operons"]:
        for ph in ("exp", "stat"):
            assert b["operons"][g][ph][2] == a["operons"][g][ph][2]
