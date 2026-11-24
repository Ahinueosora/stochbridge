"""
Module 1 - stochastic_engine (A4). Gillespie SSA of the A3 v3 CRP/lac circuit via GillesPy2 (C++ SSACSolver).

    sc = StochasticCircuit("outputs/circuit_model.yaml")            # compiles once
    traj = sc.run(InducerLevelChange(40, n_cells=1000, seed=1))     # -> io_schema.Trajectory
    traj.to_parquet("outputs/a4/traj.parquet"); traj.write_summary("outputs/a4/summary_stats.json")

Model: the reactions, parameters, initial conditions and CRP outputs are read from circuit_model.yaml (A3), the
deterministic skeleton. Inputs (TMG, glucose, crp knockout, overrides) enter through the same input functions as
A3 (reaction_network.py), so the SSA mean follows the A3 ODE wherever the ODE has one stable state.

Two modes (same species, same output columns):
    "exact" - all 8 lac reactions, internal TMG (X) simulated molecule by molecule.
    "fast"  - X held at its quasi-steady level X = k_up_eff * Y / d_x inside R1 (X relaxes in ~2 min while LacY
              changes over hours, and X is in the hundreds-thousands, so its own noise is negligible). Cuts the
              event count ~1000x for ON cells. X_count is then derived (metadata x_column = "qssa_derived").
              tests/test_ssa_reproducibility.py and outputs/a4 check fast vs exact.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))
import reaction_network as rn  # noqa: E402
from io_schema import CORE_SPECIES, Trajectory, now_utc, SCHEMA_VERSION  # noqa: E402
from perturbation import Perturbation  # noqa: E402

RATE_NAMES = ["k_rep_max", "k_unrep", "x0", "n", "k_tx_eff", "d_m", "k_tl", "d_y", "k_up_eff", "d_x", "k_pass"]


def _pow(expr: str) -> str:
    """'(a)^b' -> 'pow(a, b)' for GillesPy2 propensity strings."""
    return re.sub(r"\(([^()]*)\)\s*\^\s*([A-Za-z_][A-Za-z0-9_]*)", r"pow(\1, \2)", expr)


class StochasticCircuit:
    def __init__(self, yaml_path="outputs/circuit_model.yaml", mode="fast", include_outputs=True):
        """yaml_path: path to a circuit-model YAML, or the dict returned by Perturbation.apply() (A6 Module 2)."""
        import gillespy2 as g
        if mode not in ("fast", "exact"):
            raise ValueError("mode must be 'fast' or 'exact'")
        if isinstance(yaml_path, dict):
            self.spec = yaml_path
            raw = yaml.safe_dump(yaml_path, sort_keys=True).encode()
            self.yaml_path = Path(f"{yaml_path.get('name', 'model')}.yaml")
        else:
            self.yaml_path = Path(yaml_path)
            raw = self.yaml_path.read_bytes()
            self.spec = yaml.safe_load(raw)
        self.model_sha256 = hashlib.sha256(raw).hexdigest()
        self.genetics = dict(self.spec.get("genetics", {}))
        self.embedded = self.spec.get("perturbation")
        self.params = dict(rn.DEFAULT_PARAMS)
        self.params.update(self.spec["parameters"])
        co = self.spec.get("crp_outputs", {})
        self.outputs = co.get("targets", {}) if include_outputs else {}
        self.tpm_to_molecules = float(co.get("tpm_to_molecules", 0.0024))
        self.noise = dict(self.spec.get("noise", {}).get("parameters", {}))
        # A5: with gene-level extrinsic noise (eta2 > 0) the 21 CRP outputs are drawn outside the SSA (see _outputs)
        self.outputs_in_ssa = float(self.noise.get("output_extrinsic_cv2", 0.0)) <= 0
        self.mode = mode
        self.species = CORE_SPECIES if mode == "exact" else [s for s in CORE_SPECIES if s != "X"]

        m = g.Model(name=f"crp_lac_{mode}")
        base = self._rates(self.params, Perturbation())
        m.add_parameter([g.Parameter(name=k, expression=float(v)) for k, v in base.items()])
        sp = {s: g.Species(name=s, initial_value=0, mode="discrete") for s in self.species}
        m.add_species(list(sp.values()))
        for r in self.spec["reactions"]:
            rid, law = r["id"], r["propensity"]
            if mode == "fast":
                if rid in ("R7_uptake", "R8_tmg_loss"):
                    continue
                if rid == "R1_repress":    # internal inducer at quasi-steady state; k_pass = LacY-independent influx (IPTG)
                    law = "k_rep_max / (1 + pow((k_up_eff * Y + k_pass) / (d_x * x0), n)) * P_free"
            m.add_reaction(g.Reaction(name=rid, reactants={sp[k]: v for k, v in r["reactants"].items()},
                                      products={sp[k]: v for k, v in r["products"].items()},
                                      propensity_function=_pow(law)))
        if mode == "exact":
            m.add_reaction(g.Reaction(name="R9_passive_influx", reactants={}, products={sp["X"]: 1}, propensity_function="k_pass"))
        for gname in (self.outputs if self.outputs_in_ssa else {}):
            s = g.Species(name=f"mRNA_{gname}", initial_value=0, mode="discrete")
            m.add_species(s)
            m.add_parameter(g.Parameter(name=f"k_{gname}", expression=0.0))
            m.add_reaction([g.Reaction(name=f"tx_{gname}", reactants={}, products={s: 1}, propensity_function=f"k_{gname}"),
                            g.Reaction(name=f"deg_{gname}", reactants={s: 1}, products={}, propensity_function=f"d_m * mRNA_{gname}")])
        t0 = time.time()
        self.solver = g.SSACSolver(model=m, variable=True)
        self.compile_s = time.time() - t0
        self.model = m

    # ---------------------------------------------------------------- condition -> numbers
    def effective_params(self, pert: Perturbation):
        p = dict(self.params)
        p.update(pert.param_overrides)
        return p

    def _resolve(self, pert):
        """A Perturbation as given, or the one embedded by Perturbation.apply(); model genetics always hold."""
        if pert is None:
            if not self.embedded:
                raise ValueError("no perturbation given and none embedded in the model (use Perturbation.apply)")
            pert = Perturbation.from_dict(self.embedded)
        if self.genetics.get("crp_knockout") and not pert.crp_knockout:
            pert = pert.replace(crp_knockout=True)
        return pert

    @staticmethod
    def _rates(p, pert):
        T, G, ko = pert.TMG_uM, pert.glucose_uM, pert.crp_knockout
        return {"k_rep_max": rn.k_rep_max(p, G, ko), "k_unrep": p["k_unrep"], "x0": p["x0"], "n": p["n"],
                "k_tx_eff": p["k_tx"] * rn.c_eff(G, p, ko), "d_m": p["d_m"], "k_tl": p["k_tl"], "d_y": rn.d_y(p),
                "k_up_eff": p["k_up"] * float(rn.beta_T(T, p)) * rn.beta_G(G, p), "d_x": p["d_x"],
                "k_pass": float(p.get("k_pass", 0.0)) + p["d_x"] * 602.0 * float(getattr(pert, "IPTG_uM", 0.0))
                          * float(p.get("f_IPTG", 1.0))}      # IPTG: LacY-independent entry (A6, perturbation.InducerLevelChange)

    def condition(self, pert: Perturbation):
        """All parameter values and initial counts sent to the compiled solver for this perturbation."""
        p = self.effective_params(pert)
        rates = self._rates(p, pert)
        P_free, M, Y, X = rn.initial_state(p, pert.initial_state == "induced", pert.glucose_uM, pert.crp_knockout)
        init = {"P_free": int(round(P_free)), "P_rep": 1 - int(round(P_free)), "M": int(round(M)),
                "Y": int(round(Y)), "X": int(round(X))}
        if self.mode == "fast":
            init.pop("X")
        out_rates, out_init = {}, {}
        if self.outputs:
            tpm = rn.crp_output_tpm(pert.glucose_uM, self.outputs, p, crp_knockout=pert.crp_knockout)
            scale = p.get("output_scale", {}) or {}
            for gname, v in tpm.items():
                mol = v * self.tpm_to_molecules * float(scale.get(gname, 1.0))   # A5: medium-specific expression level
                out_rates[f"k_{gname}"] = mol * p["d_m"]
                out_init[f"mRNA_{gname}"] = int(round(mol))       # start at the steady mean of the condition
        y_on = rn.y_max(p) * rn.c_eff(pert.glucose_uM, p, pert.crp_knockout)
        y_off = y_on / rn.rho_eff(pert.glucose_uM, p, pert.crp_knockout)
        return {"rates": rates, "init": init, "out_rates": out_rates, "out_init": out_init,
                "Y_on": y_on, "Y_off": y_off, "on_threshold_Y": float(np.sqrt(y_on * y_off))}

    # ---------------------------------------------------------------- run
    def run(self, pert: Perturbation | None = None, **kw) -> Trajectory:
        """pert: a Perturbation, or None to run the one embedded by Perturbation.apply(). kw (seed, n_cells, ...) override
        run settings, e.g. sc.run(seed=rep) for replicate populations."""
        pert = self._resolve(pert)
        if kw:
            pert = pert.replace(**kw)
        cond = self.condition(pert)
        variables = {**{k: float(v) for k, v in cond["rates"].items()}, **cond["init"]}
        if self.outputs_in_ssa:
            variables.update({**{k: float(v) for k, v in cond["out_rates"].items()}, **cond["out_init"]})
        t0 = time.time()
        import gillespy2 as g
        n_pts = int(round(pert.t_end_min / pert.dt_min)) + 1
        self.solver.model.timespan(g.TimeSpan.linspace(t=pert.t_end_min, num_points=n_pts))
        res = self.solver.run(number_of_trajectories=pert.n_cells, seed=int(pert.seed), variables=variables)
        wall = time.time() - t0
        if pert.n_cells == 1 and not isinstance(res, (list, tuple)) and not hasattr(res, "__len__"):
            res = [res]
        frames = []
        tgrid = np.asarray(res[0]["time"], float)
        for i, tr in enumerate(res):
            d = {"cell_id": np.full(len(tgrid), i, dtype=np.int32), "time": tgrid}
            for s in CORE_SPECIES:
                if s == "X" and self.mode == "fast":
                    d["X_count"] = np.rint((cond["rates"]["k_up_eff"] * np.asarray(tr["Y"]) + cond["rates"]["k_pass"])
                                           / cond["rates"]["d_x"]).astype(np.int64)
                else:
                    d[f"{s}_count"] = np.rint(np.asarray(tr[s])).astype(np.int64)
            if self.outputs_in_ssa:
                for gname in self.outputs:
                    d[f"mRNA_{gname}_count"] = np.rint(np.asarray(tr[f"mRNA_{gname}"])).astype(np.int64)
            frames.append(pd.DataFrame(d))
        df = pd.concat(frames, ignore_index=True)
        if self.outputs and not self.outputs_in_ssa:
            out = self._outputs(cond, tgrid, len(res), int(pert.seed), self.effective_params(pert))
            for gname, M in out.items():
                df[f"mRNA_{gname}_count"] = M.reshape(-1)
        meta = {"schema_version": SCHEMA_VERSION,
                "units": {"time": "min", "counts": "molecules per cell"},
                "model": {"name": self.spec.get("name"), "version": self.spec.get("version"),
                          "yaml": self.yaml_path.name, "yaml_sha256": self.model_sha256},
                "perturbation": pert.to_dict(), "genetics": self.genetics, "solver": "GillesPy2 SSACSolver (direct method, C++)",
                "gillespy2_version": __import__("gillespy2").__version__,
                "mode": self.mode, "x_column": "simulated" if self.mode == "exact" else "qssa_derived",
                "crp_outputs_method": "SSA (birth-death)" if self.outputs_in_ssa else
                    "exact immigration-death sampling, per-cell gamma rate (A5 gene-level extrinsic noise)",
                "noise": self.noise,
                "seed": int(pert.seed), "rates": cond["rates"], "initial_counts": cond["init"],
                "Y_on": cond["Y_on"], "Y_off": cond["Y_off"], "on_threshold_Y": cond["on_threshold_Y"],
                "wall_time_s": wall, "created_utc": now_utc()}
        return Trajectory(df, meta)


# ---------------------------------------------------------------- A5: CRP outputs with gene-level extrinsic noise
def _outputs_exact(rates, init_mean, d_m, tgrid, n_cells, eta2, seed):
    """Each cell's mRNA_g is an immigration-death process with its own rate k_g * w_ig, w_ig ~ Gamma(mean 1, CV2 eta2),
    independent across genes and cells (PETRI-seq: no shared component). On the sampling grid the exact transition is
    m(t+dt) = Binomial(m(t), e^{-d dt}) + Poisson(k w / d * (1 - e^{-d dt})), identical in law to SSA at those times."""
    rng = np.random.default_rng(seed + 7919)
    out = {}
    dts = np.diff(tgrid)
    for gname, k in rates.items():
        w = rng.gamma(1 / eta2, eta2, n_cells) if eta2 > 0 else np.ones(n_cells)
        lam = k * w / d_m
        m = rng.poisson(lam)                                   # start at the cell's own stationary distribution
        M = np.empty((n_cells, len(tgrid)), dtype=np.int64); M[:, 0] = m
        for j, dt in enumerate(dts, 1):
            keep = np.exp(-d_m * dt)
            m = rng.binomial(m, keep) + rng.poisson(lam * (1 - keep))
            M[:, j] = m
        out[gname] = M
    return out


def _bind_outputs():
    def _outputs(self, cond, tgrid, n_cells, seed, p):
        rates = {k[2:]: v for k, v in cond["out_rates"].items()}
        return _outputs_exact(rates, None, p["d_m"], np.asarray(tgrid), n_cells,
                              float(self.noise.get("output_extrinsic_cv2", 0.0)), seed)
    StochasticCircuit._outputs = _outputs


_bind_outputs()


def run_and_save(sc: StochasticCircuit, pert: Perturbation, outdir="outputs/a4"):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    tr = sc.run(pert)
    tr.to_parquet(outdir / f"trajectory_{pert.name}.parquet")
    tr.write_summary(outdir / f"summary_stats_{pert.name}.json")
    return tr
