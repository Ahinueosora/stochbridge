"""
Module 2 - perturbation_interface (A6).

A Perturbation says what is done to the cells: environment (inducer, carbon source), genetics (gene knockouts,
parameter edits) and how the population is run (history, length, sampling, number of cells, seed).

Stable API (Student B's screen_runner, task B3, imports this module):

    from perturbation import load_perturbations, Perturbation, GeneKnockout, InducerLevelChange, CarbonSourceShift
    model = yaml.safe_load(open("outputs/a5/calibrated_circuit_model.yaml"))
    for pert in load_perturbations("config/precise1k_condition_map.yaml"):   # every defined perturbation
        modified = pert.apply(model)            # -> new circuit-model dict (input is never changed)
        sc = StochasticCircuit(modified)        # Module 1 accepts the dict; runs the embedded perturbation
        traj = sc.run(seed=rep)                 # replicate populations: change the seed

    # cheaper for many runs on one compiled engine (no recompilation):  sc = StochasticCircuit(model); sc.run(pert)

What apply() changes in the model
    parameters        param_overrides written into `parameters` (e.g. rho for operator-decoy strains, k_pass for IPTG)
    genetics          `genetics.deleted_genes`; crp / cyaA deletion -> `genetics.crp_knockout: true` (CRP-cAMP gone, the
                      A3 CRP-free basal level is used); lacI -> rho = 1 (no repression); lacY -> k_tl = 0 (no permease);
                      a CRP-output gene -> its expression set to 0
    environment       `perturbation` block: TMG / IPTG, glucose, carbon source, history, run settings
    name              "<model name>__<perturbation name>"
The A4 call style sc.run(Perturbation) is unchanged and gives the same trajectories as running the applied model.
"""
from __future__ import annotations

import copy
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import yaml

MOLECULES_PER_UM = 602.0          # 1 uM in a 1 fL E. coli cell (Milo & Phillips, Cell Biology by the Numbers)
GLUCOSE_UM_PER_G_L = 1e6 / 180.16  # 2 g/L glucose = 11,100 uM

# carbon source -> extracellular glucose seen by the circuit. The circuit's only carbon input is glucose
# (CRP activity c_crp(G) and inducer exclusion beta_G(G), A3). Non-PTS, non-glucose carbon sources therefore give G = 0
# (cAMP high: glycerol, acetate, succinate, LB/amino acids). Other PTS sugars (fructose, mannose, sucrose ...) lower
# cAMP by an amount the circuit has no independent value for, so they are refused instead of guessed.
CARBON_TO_GLUCOSE_UM = {"glucose": None, "glycerol": 0.0, "acetate": 0.0, "succinate": 0.0, "pyruvate": 0.0,
                        "LB": 0.0, "lactate": 0.0, "none": 0.0}
PTS_SUGARS_NOT_MODELLED = {"fructose", "mannose", "sucrose", "mannitol", "sorbitol", "N-acetylglucosamine", "glucosamine"}

CRP_OUTPUT_GENES = {"hflC", "tdcB", "uhpT", "cas3", "rbsC", "yahN", "fecE", "mhpF", "flhC", "flhD", "araJ", "rbsA", "nadC",
                    "galK", "rbsB", "gyrA", "rbsR", "lamB", "malK", "preA", "idnO"}   # A2-validated CRP targets (A3 outputs)

KNOCKOUT_EFFECTS = {
    "crp": "CRP-cAMP removed: c -> A3 CRP-free basal level (Kuhlman 2007, 1/240); LacI looping weakened (rho_crp_coupling)",
    "cyaA": "no cAMP -> CRP inactive: same circuit effect as crp deletion",
    "lacI": "no repressor: rho = 1, promoter always free",
    "lacY": "no permease: k_tl = 0 for LacY, no positive feedback; TMG cannot enter (IPTG still can)",
}


@dataclass
class Perturbation:
    name: str = "baseline"
    TMG_uM: float = 0.0
    glucose_uM: float = 0.0
    crp_knockout: bool = False
    param_overrides: dict = field(default_factory=dict)   # e.g. {"rho": 5} (decoy strain), {"k_pass": ...} (IPTG)
    initial_state: str = "uninduced"                       # "uninduced" (OFF history) | "induced" (ON history)
    t_end_min: float = 1200.0                              # Ozbudak protocol: 20 h
    dt_min: float = 5.0                                    # recording interval
    n_cells: int = 100
    seed: int = 1                                          # GillesPy2 needs a positive seed
    kind: str = "Perturbation"
    # A6 additions (all optional; defaults keep A4 behaviour)
    IPTG_uM: float = 0.0
    carbon_source: str = ""
    deleted_genes: list = field(default_factory=list)
    description: str = ""

    def __post_init__(self):
        if self.initial_state not in ("uninduced", "induced"):
            raise ValueError("initial_state must be 'uninduced' or 'induced'")
        if min(self.TMG_uM, self.glucose_uM, self.IPTG_uM) < 0 or self.n_cells < 1 or self.t_end_min <= 0 or self.dt_min <= 0:
            raise ValueError("bad perturbation values")

    # ------------------------------------------------------------------ serialisation
    def to_dict(self):
        return asdict(self)

    def to_json(self):
        return json.dumps(self.to_dict())

    @classmethod
    def from_dict(cls, d):
        """Rebuild from to_dict() output (subclass identity is kept in the `kind` field)."""
        return Perturbation(**{k: v for k, v in dict(d).items() if k in Perturbation.__dataclass_fields__})

    def replace(self, **kw):
        """Copy with some fields changed (e.g. seed for replicate populations)."""
        d = self.to_dict(); d.update(kw)
        return Perturbation(**d)

    # ------------------------------------------------------------------ the Module-2 contract
    def apply(self, circuit_model):
        """circuit_model: dict (parsed circuit_model.yaml) or a path to it. Returns a modified deep copy."""
        spec = yaml.safe_load(Path(circuit_model).read_text()) if isinstance(circuit_model, (str, Path)) else \
            copy.deepcopy(circuit_model)
        if "parameters" not in spec or "reactions" not in spec:
            raise ValueError("not a circuit model (needs 'parameters' and 'reactions')")
        p = spec["parameters"]
        p.update(copy.deepcopy(self.param_overrides))
        gen = spec.setdefault("genetics", {"deleted_genes": [], "crp_knockout": False})
        for g in self.deleted_genes:
            gen["deleted_genes"] = sorted(set(gen.get("deleted_genes", [])) | {g})
            if g in ("crp", "cyaA"):
                gen["crp_knockout"] = True
            elif g == "lacI":
                p["rho"] = 1.0
            elif g == "lacY":
                p["k_tl"] = 0.0
            elif g in spec.get("crp_outputs", {}).get("targets", {}):
                t = spec["crp_outputs"]["targets"][g]
                t["basal_TPM"], t["amplitude_TPM"] = 0.0, 0.0
            else:
                raise NotImplementedError(f"knockout of {g} is not in the circuit (allowed: crp, cyaA, lacI, lacY, CRP outputs)")
        if self.crp_knockout:
            gen["crp_knockout"] = True
        pert = self.to_dict()
        pert["crp_knockout"] = bool(gen.get("crp_knockout", False))
        spec["perturbation"] = pert
        spec["name"] = f"{spec.get('name', 'circuit')}__{self.name}"
        return spec


# ---------------------------------------------------------------------- subclasses
class InducerLevelChange(Perturbation):
    """Lac inducer dose. TMG (A3/A4 native): enters through LacY (positive feedback). IPTG: enters by diffusion, so the
    internal level does not depend on LacY. The engine adds a passive influx k_pass = d_x * 602 * [IPTG uM] * f_IPTG, i.e. internal
    IPTG = external concentration (602 molecules per uM in a 1 fL cell); f_IPTG = 1 assumes the same LacI affinity per molecule
    as the internal TMG pool (an assumption; at >= 20 uM IPTG the switch is saturated either way)."""

    def __init__(self, TMG_uM: float = 0.0, *, IPTG_uM: float = 0.0, f_IPTG: float | None = None, **kw):
        kw.pop("kind", None)
        po = dict(kw.pop("param_overrides", {}))
        if f_IPTG is not None:
            po["f_IPTG"] = float(f_IPTG)
        name = kw.pop("name", f"IPTG_{IPTG_uM:g}uM" if IPTG_uM and not TMG_uM else f"TMG_{TMG_uM:g}uM")
        super().__init__(name=name, TMG_uM=float(TMG_uM), IPTG_uM=float(IPTG_uM), param_overrides=po,
                         kind="InducerLevelChange", **kw)


class CarbonSourceShift(Perturbation):
    """Carbon source -> glucose seen by the circuit. CarbonSourceShift(glucose_uM) (A4 style) or
    CarbonSourceShift(carbon_source="glycerol") / ("glucose", g_per_L=2)."""

    def __init__(self, glucose_uM: float | None = None, TMG_uM: float = 0.0, *, carbon_source: str = "",
                 g_per_L: float | None = None, **kw):
        kw.pop("kind", None)
        if glucose_uM is None:
            cs = carbon_source or "none"
            if cs in PTS_SUGARS_NOT_MODELLED:
                raise NotImplementedError(f"{cs}: PTS sugar, CRP activity on it is not part of the circuit")
            if cs not in CARBON_TO_GLUCOSE_UM:
                raise NotImplementedError(f"unknown carbon source {cs}")
            glucose_uM = (g_per_L or 2.0) * GLUCOSE_UM_PER_G_L if cs == "glucose" else CARBON_TO_GLUCOSE_UM[cs]
        name = kw.pop("name", carbon_source or f"glucose_{glucose_uM:g}uM")
        super().__init__(name=name, glucose_uM=float(glucose_uM), TMG_uM=float(TMG_uM), carbon_source=carbon_source,
                         kind="CarbonSourceShift", **kw)


class GeneKnockout(Perturbation):
    """Deletion of one or more circuit genes (crp, cyaA, lacI, lacY, or a CRP-output gene)."""

    def __init__(self, gene: str | list | None = None, TMG_uM: float = 0.0, **kw):
        kw.pop("kind", None)
        kw.pop("crp_knockout", None)
        dg = kw.pop("deleted_genes", None)
        gene = gene if gene is not None else (dg or "crp")
        genes = [gene] if isinstance(gene, str) else list(gene)
        bad = [g for g in genes if g not in KNOCKOUT_EFFECTS and g not in CRP_OUTPUT_GENES]
        if bad:
            raise NotImplementedError(f"knockout of {bad} has no effect inside the circuit (allowed: crp, cyaA, lacI, lacY, "
                                      "or a CRP-output gene)")
        name = kw.pop("name", "d" + "_d".join(genes))
        crp_ko = any(g in ("crp", "cyaA") for g in genes)
        po = dict(kw.pop("param_overrides", {}))
        if "lacI" in genes:
            po["rho"] = 1.0
        if "lacY" in genes:
            po["k_tl"] = 0.0
        super().__init__(name=name, TMG_uM=float(TMG_uM), crp_knockout=crp_ko, deleted_genes=genes, param_overrides=po,
                         kind="GeneKnockout", **kw)


def OperatorDecoy(rho, TMG_uM=0.0, **kw):
    """Extra lac operators titrate LacI, lowering the repression factor (Ozbudak's rho = 5-50 strains)."""
    po = dict(kw.pop("param_overrides", {}), rho=float(rho))
    return Perturbation(name=kw.pop("name", f"decoy_rho{rho:g}"), TMG_uM=TMG_uM, param_overrides=po,
                        kind="OperatorDecoy", **kw)


KINDS = {"InducerLevelChange": InducerLevelChange, "CarbonSourceShift": CarbonSourceShift, "GeneKnockout": GeneKnockout}


def compose(*perts: Perturbation, name: str | None = None, **run_kw) -> Perturbation:
    """Several perturbations on one population (e.g. crp knockout grown on glycerol). Later ones win on conflicts;
    knockouts, overrides and inducers add up."""
    d = Perturbation().to_dict()
    for p in perts:
        q = p.to_dict()
        for k in ("TMG_uM", "IPTG_uM", "glucose_uM"):
            d[k] = max(d[k], q[k]) if k != "glucose_uM" else (q[k] if p.kind == "CarbonSourceShift" else d[k])
        d["crp_knockout"] |= q["crp_knockout"]
        d["deleted_genes"] = sorted(set(d["deleted_genes"]) | set(q["deleted_genes"]))
        d["param_overrides"].update(q["param_overrides"])
        d["carbon_source"] = q["carbon_source"] or d["carbon_source"]
    d["name"] = name or "+".join(p.name for p in perts)
    d["kind"] = "+".join(p.kind for p in perts)
    d.update(run_kw)
    return Perturbation(**d)


def _build(steps, name, run):
    parts = []
    for step in steps:
        step = dict(step)
        parts.append(KINDS[step.pop("kind")](**step))
    if len(parts) == 1:
        return parts[0].replace(name=name, **run)
    return compose(*parts, name=name, **run)


def from_spec(spec: dict, defaults: dict | None = None, control: bool = False) -> Perturbation:
    """One condition-map entry -> its Perturbation (control=True: the matched in-silico control)."""
    run = dict((defaults or {}).get("run", {}))
    run.update(spec.get("run", {}))
    steps = spec["sim_control"] if control else spec["steps"]
    p = _build(steps, spec["name"] + ("__control" if control else ""), run)
    return p.replace(description=spec.get("description", "") if not control else f"control for {spec['name']}")


def load_condition_map(config="config/precise1k_condition_map.yaml", root: Path | None = None) -> dict:
    path = Path(config)
    if not path.is_absolute() and not (root is None and path.exists()):
        path = (root or Path(__file__).parent.parent) / path
    return yaml.safe_load(path.read_text())


def load_perturbations(config="config/precise1k_condition_map.yaml", root: Path | None = None, with_controls=False):
    """Every perturbation defined in the condition map, in file order. with_controls=True returns
    [(perturbation, control, matched_dict_or_None), ...] instead."""
    cfg = load_condition_map(config, root)
    d = cfg.get("defaults", {})
    if not with_controls:
        return [from_spec(s, d) for s in cfg["perturbations"]]
    return [(from_spec(s, d), from_spec(s, d, control=True), s.get("matched")) for s in cfg["perturbations"]]
