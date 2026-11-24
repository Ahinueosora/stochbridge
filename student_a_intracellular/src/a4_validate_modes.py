"""A4 check: fast mode (quasi-steady internal TMG) vs exact SSA. Same perturbations, independent seeds.
Usage: python src/a4_validate_modes.py [n_cells]"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy import stats
sys.path.insert(0, str(Path(__file__).parent))
from ssa_simulator import StochasticCircuit
from perturbation import InducerLevelChange
from io_schema import Trajectory

OUT = Path("outputs/a4"); OUT.mkdir(parents=True, exist_ok=True)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
modes = {m: StochasticCircuit("outputs/circuit_model.yaml", mode=m, include_outputs=False) for m in ("exact", "fast")}
cases = [InducerLevelChange(40, initial_state="uninduced", n_cells=n, seed=101, name="T40_from_OFF"),
         InducerLevelChange(15, initial_state="uninduced", n_cells=n, seed=102, name="T15_from_OFF"),
         InducerLevelChange(15, initial_state="induced", n_cells=max(n // 4, 30), seed=103, name="T15_from_ON")]   # ON cells: ~30 s each exact


def first_switch(tr):
    W = tr.wide("Y"); thr = tr.meta["on_threshold_Y"]; t = W.columns.values
    st = []
    for row in W.values:
        idx = np.where(row > thr)[0] if tr.meta["perturbation"]["initial_state"] == "uninduced" else np.where(row < thr)[0]
        st.append(t[idx[0]] if len(idx) else np.inf)
    return np.array(st)


res = {"n_cells": n}
for c in cases:
    r = {}
    trs = {}
    for m, sc in modes.items():
        t0 = time.time()
        cache = OUT / f"exact_{c.name}.parquet"
        if m == "exact" and cache.exists():                      # reuse a finished exact run with the same settings
            old = Trajectory.from_parquet(cache)
            if old.meta.get("perturbation", {}).get("n_cells") == c.n_cells and old.meta.get("seed") == c.seed:
                trs[m] = tr = old
                print(c.name, "exact: reused cached run", flush=True)
            else:
                tr = sc.run(c)
        else:
            tr = sc.run(c if m == "exact" else c.replace(seed=c.seed + 1000))
        trs[m] = tr
        if m == "exact":
            tr.to_parquet(OUT / f"exact_{c.name}.parquet")
        s = tr.summary_stats()
        on_cells = tr.at(tr.times[-1]).Y_count.values
        on_cells = on_cells[on_cells > tr.meta["on_threshold_Y"]]
        r[m] = {"wall_s": time.time() - t0, "fraction_ON_20h": s["fraction_ON"],
                "mean_Y_20h": s["species"]["Y_count"]["mean"], "fano_Y_ON_cells": float(on_cells.var(ddof=1) / on_cells.mean()) if len(on_cells) > 2 else None,
                "mean_Y_ON_cells": float(on_cells.mean()) if len(on_cells) else None,
                "switched_by_20h": int(np.isfinite(first_switch(tr)).sum())}
    a, b = first_switch(trs["exact"]), first_switch(trs["fast"])
    fa, fb = np.isfinite(a).mean(), np.isfinite(b).mean()
    se = np.sqrt(fa * (1 - fa) / len(a) + fb * (1 - fb) / len(b))
    r["switched_fraction_diff"] = float(fa - fb); r["switched_fraction_diff_se"] = float(se)
    if np.isfinite(a).sum() > 5 and np.isfinite(b).sum() > 5:
        r["switch_time_KS_p"] = float(stats.ks_2samp(a[np.isfinite(a)], b[np.isfinite(b)]).pvalue)
    res[c.name] = r
    print(c.name, json.dumps(r, default=float), flush=True)
(OUT / "fast_vs_exact.json").write_text(json.dumps(res, indent=2, default=float))
