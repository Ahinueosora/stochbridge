# A3 v3: final parameterized circuit for A4

**What changed (24 Nov 2026).** The circuit is now built on network v3: A2 causal edges + RegulonDB 14.5.1 + Kuhlman 2007 + Ozbudak 2004. The wild-type model and its fit are unchanged; the knockout behaviour, the edge evidence and the SBML are new.

| | v2 | v3 |
|---|---|---|
| CRP → lac | RegulonDB + correlational | RegulonDB grade 1 (site mutation) + Kuhlman graded intervention with Δcrp control |
| crp knockout | lac output exactly 0 | CRP-free basal 1/240 of max (Kuhlman Table 2) and weaker LacI loop (ρ 170 → 22) |
| Edge labels | `a2_validated` only | `validated_by`: A2, Kuhlman2007_intervention, RegulonDB_grade1/2, Ozbudak2004_measured |
| Context | none | 85 RegulonDB regulators around the circuit (`context_edges`, not simulated) |
| SBML | GillesPy2 export, untested in COPASI | rewritten as molecule counts in minutes; checked in libSBML, GillesPy2 and COPASI |

## Outputs

| File | Content |
|---|---|
| `circuit_model.yaml` | species, reactions with propensities, parameters + sources, initial conditions (OFF and ON history), rate constants at 3 conditions, 31 circuit edges with evidence, 85 context edges, 21 CRP outputs |
| `circuit_model.xml` | the same model in SBML L3V2 at TMG 15 µM, 0 glucose, cells starting OFF (26 species, 50 reactions incl. CRP outputs). Other conditions: `build_a3_v3.write_sbml(path, p, T, G, outputs, induced, crp_knockout)` |
| `validated_edges.json` | every circuit edge and where its validation comes from |
| `hysteresis_check.png` | hysteresis (wild type + crp knockout), phase diagram, ρ estimates incl. Kuhlman, titration strains |
| `a3_checks_v3.json` | every number below |

## Checks

**Wild type (identical to v2; the new terms act only in a crp knockout):**
- Ozbudak Fig 2c: switch-on log10 RMSE 0.037 (fitted), switch-off 0.090 (held out).
- Memory window at 0 glucose: 4.1–30 µM at steady state (paper 3–30); 2.9–36 µM after 20 h.
- At 1 mM glucose: 119–871 µM (paper ~140–790).
- ON/OFF expression 170×; cusp at ρ = 9.01; ρ = 5 strain graded, ρ = 50 and 170 switch.
- Still failing: graded-response steepness at ρ = 5, Hill 2.7 vs ≈ 6. Kuhlman's looping raises it only to 3.2, so it is not the explanation.

**crp knockout (new):**
- ρ drops to 22 (CRP-assisted loop lost), but transcription falls to 1/240, so the switch-on TMG moves to ~10 mM.
- Result: LacY stays below 1 molecule per cell up to 2 mM TMG (wild type: 1,506). This matches Kuhlman (flat response in Δcrp) and Traxler 2006 (no induction without CRP).

**SBML loads and runs in both tools:**

| Check | TMG 15 µM (inside memory zone) | TMG 40 µM (above switch-on) |
|---|---|---|
| libSBML consistency | 0 errors | 0 errors |
| COPASI deterministic LacY at 20 h vs Python ODE | 10.12 vs 10.12 | 1,320.2 vs 1,320.9 |
| GillesPy2 SSA (C++), 20 cells, ON at 20 h | 4 / 20 | 11 / 20 |
| COPASI stochastic (direct method), 10 cells, LacY at 20 h | 0–519 | 1–891 |

At 15 µM the deterministic model stays OFF, but a few stochastic cells switch ON: noise-driven switching inside the memory zone, which A4 will quantify.

**Notes for A4:**
- GillesPy2 `import_SBML` makes species continuous. Set every species to `mode = "discrete"` (integer initial values) and use `SSACSolver`, or the run is silently an ODE.
- COPASI reads the file with quantity unit `#` (molecules), so its stochastic methods work directly.
- `tests/test_circuit.py`: 11 tests pass (3 new: crp knockout stays OFF, CRP-free basal = 1/240, SBML imports into GillesPy2 with 8 + 2 × 21 reactions).

## Edge evidence

- **Lac core (10 edges):** CRP → lac and LacI ⊣ lac have three independent sources each (RegulonDB grade 1, Kuhlman intervention, Ozbudak). CRP → LacI loop: Kuhlman only. Small-molecule edges (uptake, inducer → LacI, glucose → cAMP, inducer exclusion): Ozbudak measurements, Kuhlman for LacY feedback.
- **21 A2 edges:** all `a2_validated: true`. 10 also have agreeing RegulonDB grade 1–2 evidence (cas3, galK, lamB, malK, tdcB, uhpT; flhC, flhD, hflC, yahN). fecE, mhpF, nadC and gyrA disagree in sign with RegulonDB and carry only A2.
- `a2_validated` stays false for every lac edge; the brief's A2-only rule is replaced by the per-edge `validated_by` list, pending the supervisor's decision.

## Reproduce

```bash
python src/build_a3_v3.py network_v3_edges.csv outputs     # needs fitted_params.json + crp_outputs_fit.csv from run_a3.py
python -m pytest -q tests
```
Refitting from raw data still uses `run_a3.py` (needs PRECISE-1K + K-12 modulome).

---

# A3 report (v2, results unchanged in v3): CRP/lac circuit topology, rate constants and hysteresis check

**Result: the A3 model passes the Ozbudak sanity check, and it now contains all 21 A2-validated CRP edges.**

- **Phase diagram.** The model reproduces the Ozbudak et al. (2004) phase diagram:
  - switch-on thresholds, which were fitted, match within 9% (log10 RMSE 0.037);
  - switch-off thresholds, which were held out, are predicted within 23% (0.090).
- **Repression factor.** ρ ≈ 160–173 from three routes: Ozbudak's GFP (170), a free fit to his phase diagram (173), and PRECISE-1K lacY RNA in MG1655 (160). Only the RNA route is independent of Ozbudak's data.
- **The 21 A2 edges.** They are attached as CRP outputs fitted on PRECISE-1K. They predict:
  - A2's knockout effects (r = 0.85);
  - an unseen cAMP-titration experiment (r = 0.79, direction correct 93% of the time).

- **CRP → lac now has interventional support** from Kuhlman et al. 2007, a graded cAMP titration with lac ON and a Δcrp control. The edge is labelled `evidence: published_intervention`, still `a2_validated: false` (Section 4b).
- **CRP also tightens LacI repression** (Kuhlman, Ω = 12.4). The model now lowers ρ to ≈ 22 in a crp knockout. Glucose is left out of this, because letting it lower ρ breaks the phase-diagram fit (Section 4b).

Code: `src/reaction_network.py`, `src/fit_parameters.py`, `src/crp_outputs.py`, `src/run_a3.py`. Tests: `tests/test_circuit.py` (8 pass).

## 1. Why the lac edges are not A2-validated, and what would be needed

The brief asks for a network built only from A2-validated edges. A2 validated none of the lac edges, for three reasons.

**1. In the OFF state, the CRP effect on lac sits below RNA-seq's detection floor.** The lac promoter needs CRP AND no LacI. Without inducer, LacI holds the operon at 1/ρ of maximum, so CRP can only change an already tiny signal.

Predicted basal lacY is (lacY ÷ gat in induced cells) ÷ ρ × gat. Comparing that prediction with measurements:
- **Wild type:** predicted and measured match (3.7 vs 3.8 TPM, K-12 `raman` project).
- **CRP-deficient cells:** the prediction is 0.03 TPM, below RNA-seq's floor of about 1 TPM, so measured lacY stays at 1–2 TPM whatever CRP does.

The cAMP titration shows this directly (`fig_a3_lac_floor.png`):
- a cyaA knockout given 0, 0.1, 0.5 and 1 mM cAMP raises CRP activity about 40×;
- lacY stays flat at 1.1–1.6 TPM.

In A2's `crp` project even wild-type lac read only 1.6–2.5 TPM, inside the floor. A2's null result for lac is therefore a measurement limit, not evidence against the edge.

**2. Even where the effect is visible, the A2 test cannot attribute it.** In `raman`, wild type vs cyaA knockout does change basal lac (lacZ +1.08 log2, p < 0.001; lacY +0.85, p = 0.002). Conditioning on CRP activity does not remove it, because within the knockouts CRP activity varies 40× while lac does not move (it is at the floor). So the mediation test fails here too.

**3. Four of the six circuit edges are not gene → gene edges at all.** TMG inactivating LacI, LacY pumping TMG, glucose lowering cAMP, and glucose blocking LacY involve small molecules or transport. Expression-based mediation cannot test them with any amount of data. LacI → lac could be tested, but no strain varies LacI with lac intact.

**What would be needed:** RNA-seq where LacI is inactive (IPTG/TMG or ΔlacI) and CRP is perturbed (Δcrp, ΔcyaA ± cAMP), in the same experiment. Neither PRECISE-1K nor the 2,710-sample K-12 modulome has that combination.

**Best available evidence, observational:** in 19 IPTG-induced samples, log lacY tracks log gat with slope 1.07 ± 0.10 (r = 0.93). That is Ozbudak's Fig. 3b (lac output proportional to CRP) in RNA-seq, but correlational, not an intervention.

**Correction to the earlier statement "PRECISE-1K has no lac knockouts".** No sample is named as a lac mutant and none contains a lac inducer, but two backgrounds lack a functional lac operon:
- **DGF-298** (26 samples): genome-reduced; lacZYA and lacI all read zero.
- **BW25113** (160 samples): carries ΔlacZ4787 and lacI^q.

The ρ reference set now uses MG1655 only (it previously included these, inflating ρ).

## 2. Topology

| Edge | Sign | Source | A2-validated |
|---|---|---|---|
| CRP-cAMP → lac promoter | + | **Kuhlman 2007 cAMP titration (published intervention)**; RegulonDB; Ozbudak Fig. 3b; induced RNA-seq slope 1.07 | no (Section 1); published intervention (4b) |
| CRP-cAMP → LacI repression loop | + (tighter OFF) | Kuhlman 2007 Fig. 2a (published intervention) | no; in model for crp knockout only (4b) |
| LacI → lac promoter | − | RegulonDB; Ozbudak Fig. 3c | no |
| internal TMG → LacI | − | literature | not testable |
| LacY → internal TMG | + | Ozbudak Fig. 3d | not testable |
| glucose → CRP-cAMP | − | PRECISE-1K (gat 2.2× lower on glucose); Ozbudak Fig. 3a | not testable |
| glucose → LacY activity | − | inducer exclusion; Ozbudak Fig. 3e | not testable |
| CRP-cAMP → 21 targets | ± | A2 causal-mediation test | **yes** (Section 5) |

Core species: promoter free/repressed, lac mRNA, LacY, internal TMG. Eight reactions (R1–R8), plus one mRNA species per CRP output with its own transcription and decay reactions. At steady state the lac core reduces exactly to Ozbudak's u = a(1 + u²)/(ρ + u²).

## 3. What the data constrain

### 3.1 ρ from RNA-seq

ρ = (lac ÷ gat) in fully IPTG-induced cells ÷ (lac ÷ gat) in uninduced cells. gat is the CRP-only promoter that served as Ozbudak's red reporter, so dividing by it cancels CRP differences between samples. Induced groups with raw lacY from 100 to 1,500 TPM all give lacY ÷ gat between 0.12 and 0.19.

| Gene | Uninduced reference | ρ | 95% CI |
|---|---|---|---|
| lacY | PRECISE-1K MG1655 (n = 560) | 160 | 144–288 |
| lacY | same batch, K-12 srna_oe WT (n = 3) | 66 | 32–118 |
| lacZ | PRECISE-1K MG1655 | 435 | 403–724 |
| lacZ | same batch | 92 | 69–189 |

Every estimate is far above 9, the minimum for a switch, so the data predict bistability. lacY is the more trustworthy gene: several induced strains carry lac-derived plasmid sequence that may add lacZ reads. The same-batch reference has only 3 samples sitting at the floor.

### 3.2 Which ρ the phase diagram picks

A free fit of ρ to all 14 Fig. 2c thresholds gives **173**. Fixing ρ at the other values makes the held-out switch-off fit worse:

| ρ used | Switch-off error (log10 RMSE) |
|---|---|
| 170 (model) | 0.090 |
| 66 | 0.38 |
| 435 | 0.30 |

### 3.3 CRP and glucose (PRECISE-1K, MG1655)

- CRP activity (gat) is 2.2× lower in glucose-grown cells (338) than in other carbon sources (227), p = 10⁻³². Use this for the A6 carbon-source perturbation.
- Basal lacY is only 1.5× lower, with a slope against gat of 0.30 (Ozbudak expects 1). This compression is the floor effect from Section 1.

## 4. Lac core parameters and hysteresis check

| Parameter | Value | Source |
|---|---|---|
| ρ | 170 | Section 3 |
| inducer–LacI cooperativity n | 2 | assumption (gives cusp at ρ = 9) |
| max / basal LacY | 1,506 / 9 molecules | set: transcription 0.58/min, mRNA half-life 1.5 min, burst 20, doubling 90 min |
| LacI off-rate | 0.1/min | set; affects noise only (A4) |
| uptake | k_up · TMG^0.6, k_up = 1.84 | exponent from Ozbudak Fig. 3d; scale fitted (switch-on = 30 µM) |
| CRP at saturating glucose / half-effect | 0.20 / 3.9 µM | Ozbudak Fig. 3e / fitted to Fig. 2c switch-on |
| LacY activity at saturating glucose / half-effect | 0.65 / 7.4 µM | Ozbudak Fig. 3e / fitted to Fig. 2c switch-on |

| Check | Paper | Model | Pass |
|---|---|---|---|
| Bistable window, 0 glucose (steady state) | 3–30 µM | 4.1–30 µM | yes |
| Same, 20 h protocol | 3–30 µM | 2.9–36 µM | yes |
| ON ÷ OFF expression | >100× | 141× | yes |
| Window at 1 mM glucose | ~140–790 µM | 119–871 µM | yes |
| Cusp | ρ ≈ 9 | 9.01 | yes |
| ρ = 50 strain / ρ = 5 strain | switch / graded | switch / graded | yes |
| Steepness of graded response (Hill) | ≈ 6 (paper's curve); the 5 points allow 2.7–6.0 | 2.7 log-grid / 2.3 same fit as paper (3.2 with looping; unchanged with titration) | **borderline** (Section 4c) |
| Half-induction of the ρ = 5 strain (1 mM glucose) | ~160 µM | ~5 µM | **no** (Section 4c) |

The 20 h simulation is wider than the steady state because cells near a threshold switch slowly. With 200 h it converges (30.2 / 4.07 µM).

## 4b. Evidence from Kuhlman 2007, and CRP-dependent ρ

Full analysis: `kuhlman2007/kuhlman_report.md`. The data were read exactly from the vector figures in the paper's PDF.

**CRP → lac, tested by intervention.** Kuhlman's cells cannot make or break cAMP (ΔcyaA ΔcpdA ΔlacY), so added cAMP sets CRP activity, and 1 mM IPTG keeps LacI off.
- lac rises 201× from 0 to 1 mM cAMP.
- In Δcrp, the response is flat (slope −0.014 ± 0.006).
- Output is close to proportional to CRP occupancy (Hill 1.15, CI 1.03–1.28; paper 1.0 ± 0.1), which is what reaction R3 assumes.

**CRP strengthens LacI's repression loop.** Ω = 12.4 (CI 10.2–15.4), so ρ − 1 grows about 10× from no CRP to full CRP. In the model:
ρ_eff = 1 + (ρ − 1)(1 + L(θ)) / (1 + L(θ_ref)), with L(θ) = L0(1 + (Ω − 1)θ), L0 = 3.91 and θ_ref = 0.76.

| Where ρ follows CRP (`rho_crp_coupling`) | ρ at 1 mM glucose | ρ, crp knockout | Switch-on error (fit) | Switch-off error (held out) |
|---|---|---|---|---|
| nowhere (old model) | 170 | 170 | 0.037 | 0.090 |
| **crp knockout only (default)** | 170 | **22** | **0.037** | **0.090** |
| glucose too, Kuhlman strength | 52 | 22 | 0.58 | 0.32 |
| glucose too, weak (θ_ref = 0.05) | 128 | 117 | 0.11 | 0.16 |

**Why glucose is left out.** If glucose lowers ρ, the switch-on threshold stops rising with glucose. Lower CRP cuts both the drive and ρ, and the two effects cancel, so the model can no longer fit Ozbudak Fig. 2c. Every step up in coupling makes the fit worse. Kuhlman's own wild-type data agree with leaving it out: MG1655 on glucose keeps an IPTG fold change of 1,200, close to the ~1,500 seen with saturating cAMP. Ordinary glucose barely weakens the loop; removing CRP genetically does.

One hint points the other way. Ozbudak's ON/OFF threshold ratio shrinks from 10 (no glucose) to 5.6 (1 mM). With fixed ρ the model gives 7.3 at all glucose. That gap is within the ±15% digitizing error on each threshold.

**Consequence for A6.** A crp-knockout perturbation should use ρ ≈ 22 and c = 0. One mismatch remains: R3 makes lac transcription exactly 0 when c = 0, but Kuhlman measures about 1/200 of maximum without CRP. Add a CRP-independent floor to R3 if A6 needs lac levels in Δcrp.

**Looping does not explain the steepness gap.** Replacing the promoter with Kuhlman's looping form raises the ρ = 5 graded Hill coefficient from 2.7 to only 3.2; Ozbudak measured ≈ 6.

## 4c. Ozbudak's low-ρ strains: titration and steepness

Full analysis: `steepness/steepness_report.md`. Ozbudak lowered ρ by titration: plasmids with 4 or 25 copies of a lac promoter soak up LacI.

- **Titration modelled explicitly.** Its two constants are fixed exactly by ρ = 50 and ρ = 5. The plasmid operators come out binding LacI 22.5× more weakly than the chromosomal promoter, so they act as a linear sponge. Steepness is unchanged (Hill 2.29 vs 2.30). Titration does not explain the gap.
- **The gap is smaller than reported.** Fig. 4e has only 5 points; they allow Hill 2.7–6.0 (best 4.0). The model gives 2.3 with the same fit. "≈ 6" is the paper's fitted curve.
- **New mismatch: where the ρ = 5 strain switches on.** The model's half-induction is at ~5 µM TMG; the data's is ~160 µM. Thresholds in the plasmid strains fall much less than the model predicts. This affects A6's OperatorDecoy validation.

## 5. The 21 A2 edges as CRP outputs

Each target's mRNA is a Hill function of the same CRP activity that drives the lac promoter:
- **activated targets:** TPM = basal + amplitude × c^h / (K^h + c^h);
- **repressed targets:** the same with 1 − Hill.

Direction comes from A2. c is gat TPM divided by the median of MG1655 non-glucose samples (2,897 TPM), so c = 1 is the no-glucose state of the lac model and a crp knockout is c = 0. Fitted on 1,009 PRECISE-1K samples (DGF-298 excluded). Full table: `crp_outputs_fit.csv`.

| Check | Result |
|---|---|
| Predicted vs A2-measured knockout effect (21 genes) | r = 0.85, p = 10⁻⁶; predicted size ≈ 0.5× measured |
| Out of sample: `raman` cAMP titration, 21 genes × 4 doses vs WT | r = 0.79, p = 10⁻¹⁸; direction correct 93% |
| Median R² across all 1,009 samples | 0.22 (other regulators, e.g. MalT, RbsR, FlhDC inputs, act too) |

Notes:
- The model under-predicts knockout effects by about half. Fitting across all conditions dilutes the CRP dependence, because CRP is only one of each target's regulators.
- **Weak fits (R² < 0.06): araJ, fecE, mhpF, rbsR.** mhpF fits flat: across conditions it rises with CRP, but A2 shows the knockout raises it. That matches A2's flag that mhpF, fecE, nadC and gyrA are causal but likely indirect. They are kept, with `fit_quality: weak` in the YAML.
- The outputs do not feed back on lac. They are readouts for A6's DiMA comparisons and for Student B.

## 6. Limitations

- Lac edges rest on prior knowledge (Section 1). This is a deliberate deviation from the brief, for documented reasons.
- The Fig. 2c thresholds were digitized by eye (±15%). The glucose-response shapes are assumed hyperbolic.
- The induced RNA-seq samples are from lacI^q/plasmid strains in another project.
- Units come from chosen expression rates (about 1,500 LacY when fully induced). Outputs convert TPM to molecules with 0.0024 per TPM (about 2,400 mRNAs per cell, an assumption). A5 should recalibrate both.

## 7. Hand-off to A4

- `circuit_model.yaml`: full spec, including `crp_outputs`.
- `circuit_model_TMG15_G0.xml`: SBML at 15 µM TMG, 0 glucose, with the 21 outputs.
- `reaction_network.to_gillespy2(p, T, G, outputs=...)`: a GillesPy2 smoke test runs.

Next in A4: run ~1,000 cells from each history, check for two-peaked LacY histograms inside 4–30 µM, and test the graded-response steepness mismatch.
