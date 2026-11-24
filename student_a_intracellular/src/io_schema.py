"""
Module 1 output contract: the Trajectory format (StochBridge, Student A -> Student B Module 4).

A Trajectory is ONE run of ONE perturbation: many cells, many time points, long format (one row per cell per time).

Columns (in this order):
    cell_id            int32    0 .. n_cells-1 (one stochastic cell = one agent)
    time               float64  minutes since the perturbation started (same grid for every cell)
    P_free_count       int64    lac promoter free of LacI (0/1)
    P_rep_count        int64    lac promoter bound by LacI (0/1)
    M_count            int64    lacZYA mRNA molecules
    Y_count            int64    LacY permease molecules (the lac switch read-out)
    X_count            int64    internal TMG molecules (see metadata x_column: "simulated" or "qssa_derived")
    mRNA_<gene>_count  int64    optional, one per A2-validated CRP target (21 genes), only if include_outputs

File-level metadata (parquet key-value, all JSON strings):
    schema_version, units, model (name, version, yaml sha256), perturbation (full JSON), solver, mode, x_column,
    seed, created_utc

Change policy: core columns and their order are frozen at schema_version 1.x. Adding optional columns is a minor
version bump; renaming/removing/retyping anything is a major bump and must be flagged to Student B first.
See shared/schemas/trajectory_schema.md.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

SCHEMA_VERSION = "1.0.0"
CORE_SPECIES = ["P_free", "P_rep", "M", "Y", "X"]
CORE_COLUMNS = ["cell_id", "time"] + [f"{s}_count" for s in CORE_SPECIES]
UNITS = {"time": "min", "*_count": "molecules per cell"}


@dataclass
class Trajectory:
    data: pd.DataFrame                       # long format, columns as above
    meta: dict = field(default_factory=dict)

    # ---------------------------------------------------------------- views
    @property
    def species_columns(self):
        return [c for c in self.data.columns if c.endswith("_count")]

    @property
    def n_cells(self):
        return int(self.data.cell_id.nunique())

    @property
    def times(self):
        return np.sort(self.data.time.unique())

    def at(self, t):
        """All cells at the recorded time closest to t (one row per cell)."""
        tt = self.times[np.argmin(np.abs(self.times - t))]
        return self.data[self.data.time == tt]

    def wide(self, species="Y"):
        """cells x time matrix of one species."""
        return self.data.pivot(index="cell_id", columns="time", values=f"{species}_count")

    # ---------------------------------------------------------------- stats
    def summary_stats(self, t=None, on_threshold=None):
        """mean, variance, Fano factor per species across cells at time t (default: last time point),
        plus the same averaged over the second half of the run, and the ON fraction for LacY."""
        snap = self.at(self.times[-1] if t is None else t)
        late = self.data[self.data.time >= self.times[-1] / 2]
        out = {"time_min": float(snap.time.iloc[0]), "n_cells": self.n_cells, "species": {}}
        for c in self.species_columns:
            v = snap[c].to_numpy(float)
            vl = late[c].to_numpy(float)
            out["species"][c] = {"mean": float(v.mean()), "variance": float(v.var(ddof=1)),
                                 "fano": float(v.var(ddof=1) / v.mean()) if v.mean() > 0 else None,
                                 "second_half_mean": float(vl.mean()), "second_half_variance": float(vl.var(ddof=1)),
                                 "second_half_fano": float(vl.var(ddof=1) / vl.mean()) if vl.mean() > 0 else None}
        if on_threshold is None:
            on_threshold = self.meta.get("on_threshold_Y")
        if on_threshold is not None:
            out["Y_on_threshold"] = float(on_threshold)
            out["fraction_ON"] = float((snap.Y_count > on_threshold).mean())
        out["bimodality_Y"] = bimodality_of_counts(snap.Y_count.to_numpy())
        return out

    # ---------------------------------------------------------------- io
    def to_parquet(self, path):
        import pyarrow as pa
        import pyarrow.parquet as pq
        validate(self.data)
        table = pa.Table.from_pandas(self.data, preserve_index=False)
        md = {k.encode(): json.dumps(v, default=float).encode() for k, v in self.meta.items()}
        table = table.replace_schema_metadata({**(table.schema.metadata or {}), **md})
        pq.write_table(table, str(path), compression="zstd")
        return Path(path)

    @classmethod
    def from_parquet(cls, path):
        import pyarrow.parquet as pq
        table = pq.read_table(str(path))
        md = {k.decode(): json.loads(v.decode()) for k, v in (table.schema.metadata or {}).items()
              if not k.decode().startswith("pandas")}
        df = table.to_pandas()
        validate(df)
        return cls(df, md)

    def write_summary(self, path, **kw):
        s = self.summary_stats(**kw)
        s["perturbation"] = self.meta.get("perturbation")
        Path(path).write_text(json.dumps(s, indent=2, default=float))
        return s


def validate(df: pd.DataFrame):
    """Raise if a DataFrame breaks the Trajectory contract."""
    cols = list(df.columns)
    if cols[: len(CORE_COLUMNS)] != CORE_COLUMNS:
        raise ValueError(f"core columns must be {CORE_COLUMNS} in this order, got {cols[:len(CORE_COLUMNS)]}")
    extra = cols[len(CORE_COLUMNS):]
    bad = [c for c in extra if not (c.startswith("mRNA_") and c.endswith("_count"))]
    if bad:
        raise ValueError(f"unknown columns {bad}; optional columns must be mRNA_<gene>_count")
    if str(df.cell_id.dtype) != "int32" or str(df.time.dtype) != "float64":
        raise ValueError("cell_id must be int32 and time float64")
    for c in cols[2:]:
        if str(df[c].dtype) != "int64":
            raise ValueError(f"{c} must be int64")
        if (df[c] < 0).any():
            raise ValueError(f"{c} has negative counts")
    return True


def bimodality(x):
    """Sarle's bimodality coefficient (> 0.555 suggests two peaks) and Hartigan's dip test on log10(Y + 1)."""
    from scipy import stats
    x = np.asarray(x, float)
    n = len(x)
    if n < 4 or np.std(x) == 0:
        return {"sarle_bc": None, "dip": None, "dip_p": None}
    g, k = stats.skew(x), stats.kurtosis(x)
    bc = (g**2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    try:
        import diptest
        dip, p = diptest.diptest(x)
    except Exception:  # noqa: BLE001
        dip, p = None, None
    return {"sarle_bc": float(bc), "dip": None if dip is None else float(dip), "dip_p": None if p is None else float(p)}


def bimodality_of_counts(y, seed=0):
    """Bimodality of log10(LacY + 1) after spreading each integer count uniformly over [y, y + 1)
    (de-quantisation). Without it the dip test reads the spikes at 0, 1, 2 ... molecules as extra peaks."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y, float)
    return bimodality(np.log10(y + rng.uniform(0, 1, len(y)) + 1))


def now_utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")
