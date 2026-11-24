"""Writes notebooks/03_vilar_reproduction.ipynb (run from student_a_intracellular/)."""
import nbformat as nbf
nb = nbf.v4.new_notebook()
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
nb.cells = [
md("""# A4: all-or-none lac induction with our own circuit (Vilar, Guet & Leibler 2003)

**Question.** Does our Gillespie simulation of the A3 v3 circuit give *all-or-none* induction (two separate groups of cells,
each cell switching suddenly at its own time), and not just one noisy group?

**Vilar et al. 2003 (J Cell Biol 161:471), Fig. 3:** (a) at a sub-optimal inducer dose single cells switch abruptly at
random times while the average of 2,000 cells rises smoothly, and a population curve with slower-growing induced cells
lies below it; (b) at a maintenance dose neither induced nor uninduced cells switch; (c) at a high dose every cell switches.

**Here:** our `circuit_model.yaml` (A3 v3), not Vilar's parameters. Doses are scaled to our memory zone (4.1–30 µM TMG).
Simulations: `src/a4_vilar.py` (cached in `outputs/a4/`); this notebook only reads and plots them.
To rerun everything (~20 min): uncomment the `simulate()` line."""),
code("""import sys, json
from pathlib import Path
ROOT = Path.cwd() if (Path.cwd() / "src").exists() else Path.cwd().parent
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
import a4_vilar as V
# V.simulate()        # uncomment to rerun all simulations
R = json.loads((V.OUT / "vilar_results.json").read_text())
plt = V._style()
%matplotlib inline"""),
md("## 1. Memory survives noise: SSA hysteresis (1,000 cells per point, both histories)"),
code("""V.fig_hysteresis(plt, save=False)
pd.DataFrame([R["ssa_hysteresis"]]).T.rename(columns={0: "value"})"""),
md("## 2. Vilar Fig. 3: cells switch suddenly, at different times (cellular + population levels)"),
code("""V.fig_vilar3(plt, save=False)
pd.DataFrame({k: {kk: v[kk] for kk in ["TMG_uM", "history", "n_cells", "fraction_switched_30h", "switch_time_median_min",
                                        "final_cell_avg_Y", "final_population_avg_Y"]} for k, v in R["vilar"].items()}).T"""),
code("""a = R["vilar"]["vilar3a_sub"]
print(f"one cell goes from 10% to 90% ON in {a['single_cell_rise_10_90_median_min']:.0f} min (median);")
print(f"the 2,000-cell average takes {a['population_rise_10_90_min']:.0f} min: the average is smooth only because cells switch at different times.")"""),
md("## 3. Bimodal, not just noisy: wild type vs a graded (ρ = 5) strain with the same average"),
code("""V.fig_bimodality(plt, save=False)
rows = []
for t, b in R["bimodality_wt"].items():
    rows.append({"hours": int(t) / 60, "wild type dip p": b["dip_p"], "wild type Sarle BC": b["sarle_bc"], "wild type % ON": 100 * b["fraction_ON"],
                 "rho=5 dip p": R["rho5_control"]["bimodality_rho5"][t]["dip_p"], "rho=5 Sarle BC": R["rho5_control"]["bimodality_rho5"][t]["sarle_bc"]})
pd.DataFrame(rows).round(3)"""),
md("## 4. Molecular level: one cell switching (exact simulation, every reaction simulated)"),
code("""V.fig_molecular(plt, save=False);"""),
md("## 5. Checks: fast mode vs exact SSA, and the noise knob k_unrep"),
code("""fe = json.loads((V.OUT / "fast_vs_exact.json").read_text())
pd.DataFrame({k: {"exact % ON": 100 * v["exact"]["fraction_ON_20h"], "fast % ON": 100 * v["fast"]["fraction_ON_20h"],
                  "difference (pp)": 100 * v["switched_fraction_diff"], "± 1 SE (pp)": 100 * v["switched_fraction_diff_se"],
                  "switch-time KS p": v.get("switch_time_KS_p"), "exact ON-cell mean LacY": v["exact"]["mean_Y_ON_cells"],
                  "fast ON-cell mean LacY": v["fast"]["mean_Y_ON_cells"]} for k, v in fe.items() if k != "n_cells"}).T.round(3)"""),
code("""pd.DataFrame(R["k_unrep_sensitivity"]).pivot_table(index=["TMG_uM", "history"], columns="k_unrep", values="fraction_ON_20h").round(3)"""),
md("## 6. Conclusion\n\nSee `outputs/a4/a4_report.md` for the numbers in words."),
]
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "notebooks/03_vilar_reproduction.ipynb")
print("ok")
