"""
A3 - CRP/lac circuit: species, reactions, propensities (StochBridge, Student A).

Topology (edge provenance in config/crp_lac_circuit.yaml):
    TMG_ext --LacY--> TMG_in --| LacI --| P_lac --> lacY mRNA --> LacY   (positive feedback)
    glucose --| CRP-cAMP --> P_lac                                       (CRP route)
    glucose --| LacY activity                                            (inducer exclusion)

Species (molecule counts per cell):
    P_free  promoter not bound by LacI (0 or 1)
    P_rep   promoter bound by LacI (0 or 1)          P_free + P_rep = 1
    M       lacZYA mRNA
    Y       LacY permease
    X       internal TMG

Reactions:
    R1  P_free -> P_rep        k_rep_max / (1 + (X/x0)^n) * P_free      LacI binds (inducer-free LacI)
    R2  P_rep  -> P_free       k_unrep * P_rep                          LacI falls off
    R3  P_free -> P_free + M   k_tx * c_crp(G) * P_free                 transcription, scaled by CRP-cAMP
    R4  M      -> 0            d_m * M                                  mRNA decay
    R5  M      -> M + Y        k_tl * M                                 translation
    R6  Y      -> 0            d_y * Y                                  dilution by growth
    R7  0      -> X            k_up * beta_T(T) * beta_G(G) * Y         TMG uptake through LacY
    R8  X      -> 0            d_x * X                                  TMG efflux + dilution

Repression factor: rho = 1 + k_rep_max / k_unrep (max / basal expression).
CRP helps LacI close its repression loop (Kuhlman et al. 2007 PNAS, Omega ~ 12), so rho depends on CRP:
    rho(G) = 1 + (rho - 1) * (1 + L(theta(G))) / (1 + L(theta_ref)),   L(theta) = L0 (1 + (Omega - 1) theta)
    theta(G) = crp_occ_ref * c_crp(G)   (CRP occupancy; c_crp = 1 is the no-glucose state)
p["rho"] is the value with normal CRP. rho_crp_coupling picks where this applies: "knockout" (default, CRP
removed only), "glucose" (also glucose), or False (fixed rho).
At steady state the model collapses exactly to Ozbudak's form
    u = a (1 + u^n) / (rho(G) + u^n),  u = X/x0,  a = Y_max c_crp k_up beta_T beta_G / (d_x x0)
"""
from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------- parameters
# Time unit: minutes. Values marked "fit" are overwritten by fit_parameters.py.
DEFAULT_PARAMS = {
    # promoter (rho = 1 + k_rep_max / k_unrep)
    "rho": 170.0,            # fit/prior: repression factor at no glucose (c_crp = 1)
    # CRP-assisted LacI looping (Kuhlman 2007 Fig 2a, refit to vector-extracted points; kuhlman2007/kuhlman_report.md)
    "rho_crp_coupling": "knockout",  # "knockout" | "glucose" | False (see rho_eff)
    "loop_L0": 3.91,         # looping strength with no CRP
    "loop_omega": 12.4,      # CRP-loop cooperativity (paper 10.3; refit 12.4, CI 10.2-15.4)
    "crp_occ_ref": 0.76,
    # lac output with no CRP at all, as a fraction of full CRP activity (Kuhlman 2007 Table 2: 1/f_cAMP = 1/240)
    "crp_free_basal": 1 / 240,     # CRP occupancy at c_crp = 1: wt on glycerol ~ dcyaA dcpdA + 1 mM cAMP (1000/(1000+320))
    "k_unrep": 0.1,          # 1/min, complete LacI dissociation (A4 sensitivity knob)
    "n": 2.0,                # inducer-LacI cooperativity (model assumption)
    "x0": 500.0,             # molecules, internal TMG half-inactivating LacI (sets unit of X)
    # expression
    "k_tx": 0.58,            # 1/min, max transcription at full CRP activity
    "d_m": 0.46,             # 1/min, lac mRNA half-life ~1.5 min
    "k_tl": 9.2,             # 1/min per mRNA (burst size k_tl/d_m = 20)
    "doubling_min": 90.0,    # min, growth in succinate minimal medium
    # transport
    "d_x": 0.5,              # 1/min, internal TMG loss
    "k_up": None,            # fit: uptake scale (per LacY per min)
    "uptake_form": "power",  # 'power' (default) or 'hill'
    "uptake_h": 0.6,         # Ozbudak Fig 3d power-law exponent
    "uptake_K": 680.0,       # uM, Ozbudak Fig 3d half-saturation
    # glucose responses (endpoints from Ozbudak Fig 3e; half-points fit to Fig 2c)
    "crp_min": 0.20,         # CRP activity at saturating glucose (80% reduction)
    "crp_K": None,           # fit, uM glucose
    "excl_min": 0.65,        # LacY activity at saturating glucose (35% reduction)
    "excl_K": None,          # fit, uM glucose
}


def d_y(p):
    return np.log(2) / p["doubling_min"]


def y_max(p):
    """Steady LacY with promoter always free and full CRP activity."""
    return p["k_tx"] * p["k_tl"] / (p["d_m"] * d_y(p))


def k_rep_max(p, G=0.0, crp_knockout=False):
    return (rho_eff(G, p, crp_knockout) - 1.0) * p["k_unrep"]


# ------------------------------------------------------- input functions
def c_crp(G, p):
    """CRP-cAMP activity (0-1) vs extracellular glucose G (uM)."""
    return p["crp_min"] + (1 - p["crp_min"]) / (1 + G / p["crp_K"])


def c_eff(G, p, crp_knockout=False):
    """CRP activity seen by the lac promoter. Wild type: c_crp(G). crp/cyaA knockout: the CRP-free
    basal level (Kuhlman 2007: delta-crp keeps ~1/240 of full activity), not exactly 0."""
    return p.get("crp_free_basal", DEFAULT_PARAMS["crp_free_basal"]) if crp_knockout else c_crp(G, p)


def rho_eff(G, p, crp_knockout=False):
    """Repression factor: CRP strengthens LacI looping, so less CRP -> smaller rho.
    rho_crp_coupling: "knockout" (default) applies it only when CRP is removed (crp/cyaA knockout);
    "glucose" also lowers rho as glucose lowers c_crp; False keeps rho fixed.
    "knockout" is the default because the glucose form cannot reproduce Ozbudak Fig 2c (switch-on RMSE
    0.037 -> 0.58), and wild type on glucose keeps a near-full ON/OFF ratio (Kuhlman Table 1: 1,200 vs ~1,500)."""
    mode = p.get("rho_crp_coupling", False)
    if not mode or (mode == "knockout" and not crp_knockout):
        return p["rho"]
    L = lambda theta: p["loop_L0"] * (1 + (p["loop_omega"] - 1) * theta)
    theta = 0.0 if crp_knockout else p["crp_occ_ref"] * c_crp(G, p)
    return 1 + (p["rho"] - 1) * (1 + L(theta)) / (1 + L(p["crp_occ_ref"]))


def beta_G(G, p):
    """Fraction of LacY active vs glucose (inducer exclusion)."""
    return p["excl_min"] + (1 - p["excl_min"]) / (1 + G / p["excl_K"])


def beta_T(T, p):
    """Uptake per active LacY vs extracellular TMG T (uM).
    'power': (T/1 uM)^0.6, the Ozbudak Fig 3d power-law fit over the threshold range (default).
    'hill' : T^h/(K^h+T^h) with K = 680 uM (adds saturation; tested as an alternative)."""
    h, K = p["uptake_h"], p["uptake_K"]
    T = np.asarray(T, float)
    if p.get("uptake_form", "power") == "power":
        return T**h
    return T**h / (K**h + T**h)


def drive(T, G, p, crp_knockout=False):
    """Dimensionless drive a: the only combination of inputs the steady state sees."""
    return y_max(p) * c_eff(G, p, crp_knockout) * p["k_up"] * beta_T(T, p) * beta_G(G, p) / (p["d_x"] * p["x0"])


# ------------------------------------------------------- deterministic model
def rhs(t, y, T, G, p, crp_knockout=False):
    P_free, M, Y, X = y
    P_rep = 1.0 - P_free
    u = max(X, 0.0) / p["x0"]
    r_rep = k_rep_max(p, G, crp_knockout) / (1 + u ** p["n"]) * P_free
    r_unrep = p["k_unrep"] * P_rep
    dP = r_unrep - r_rep
    dM = p["k_tx"] * c_eff(G, p, crp_knockout) * P_free - p["d_m"] * M
    dY = p["k_tl"] * M - d_y(p) * Y
    dX = p["k_up"] * beta_T(T, p) * beta_G(G, p) * Y - p["d_x"] * X
    return [dP, dM, dY, dX]


def initial_state(p, induced: bool, G=0.0, crp_knockout=False):
    """Uninduced: repressed promoter, basal LacY. Induced: fully induced cell."""
    c = c_eff(G, p, crp_knockout)
    if induced:
        Y = y_max(p) * c
        return [1.0, p["k_tx"] * c / p["d_m"], Y, 50 * p["x0"]]       # steady mRNA (also valid when k_tl = 0, lacY KO)
    rho = rho_eff(G, p, crp_knockout)
    Y = y_max(p) * c / rho
    return [1.0 / rho, p["k_tx"] * c / (rho * p["d_m"]), Y, 0.0]


def simulate(T, G, p, induced: bool, t_end=20 * 60, crp_knockout=False):
    """Grow a cell 20 h (Ozbudak protocol) from a defined history; return final LacY."""
    from scipy.integrate import solve_ivp
    sol = solve_ivp(rhs, (0, t_end), initial_state(p, induced, G, crp_knockout), args=(T, G, p, crp_knockout),
                    method="LSODA", rtol=1e-7, atol=1e-9)
    return sol.y[2, -1]


# ------------------------------------------------------- analytic steady states
def a_of_u(u, rho, n=2.0):
    """Drive that holds the cell at internal-inducer level u (steady state solved for a)."""
    return u * (rho + u**n) / (1 + u**n)


def switching_drives(rho, n=2.0):
    """(a_off, a_on): saddle-node drives. None if monostable (rho below the cusp)."""
    u = np.logspace(-3, 3, 20000)
    a = a_of_u(u, rho, n)
    da = np.diff(a)
    turn = np.where(np.sign(da[:-1]) != np.sign(da[1:]))[0]
    if len(turn) < 2:
        return None
    a_on, a_off = a[turn[0] + 1], a[turn[1] + 1]   # local max, local min
    return a_off, a_on


def tmg_for_drive(a_target, G, p, crp_knockout=False):
    """Invert drive(T) for T at fixed glucose (drive is monotone in T)."""
    h, K = p["uptake_h"], p["uptake_K"]
    if p.get("uptake_form", "power") == "power":
        return (a_target / drive(1.0, G, p, crp_knockout)) ** (1 / h)
    scale = drive(1e12, G, p, crp_knockout)          # saturated uptake
    f = a_target / scale
    if f >= 1:
        return np.inf
    return K * (f / (1 - f)) ** (1 / h)


def thresholds_tmg(G, p, crp_knockout=False):
    sd = switching_drives(rho_eff(G, p, crp_knockout), p["n"])
    if sd is None:
        return None
    return tmg_for_drive(sd[0], G, p, crp_knockout), tmg_for_drive(sd[1], G, p, crp_knockout)   # (switch-off, switch-on)


# ------------------------------------------------------- A2 CRP outputs
def crp_output_tpm(G, outputs: dict, p, crp_knockout=False):
    """Steady-state expression (TPM) of each A2-validated CRP target at glucose G.
    outputs: {gene: dict(basal_TPM, amplitude_TPM, K_crp, hill_n, direction)}.
    CRP activity is the same c_crp(G) that drives the lac promoter (0 for a crp knockout)."""
    c = 0.0 if crp_knockout else c_crp(G, p)
    out = {}
    for gname, o in outputs.items():
        cc = max(c, 1e-6)
        act = cc ** o["hill_n"] / (o["K_crp"] ** o["hill_n"] + cc ** o["hill_n"])
        f = act if o["direction"] == "activated" else 1 - act
        out[gname] = o["basal_TPM"] + o["amplitude_TPM"] * f
    return out


# ------------------------------------------------------- stochastic export
def to_gillespy2(p, T, G, name="crp_lac", outputs: dict | None = None, tpm_to_molecules=0.0024):
    """GillesPy2 model for A4 (custom propensities; inputs T, G fixed per run)."""
    import gillespy2 as g
    m = g.Model(name=name)
    par = {
        "k_rep_max": k_rep_max(p, G), "k_unrep": p["k_unrep"], "x0": p["x0"], "n": p["n"],
        "k_tx_eff": p["k_tx"] * c_crp(G, p), "d_m": p["d_m"], "k_tl": p["k_tl"], "d_y": d_y(p),
        "k_up_eff": p["k_up"] * float(beta_T(T, p)) * beta_G(G, p), "d_x": p["d_x"],
    }
    m.add_parameter([g.Parameter(name=k, expression=float(v)) for k, v in par.items()])
    s = {k: g.Species(name=k, initial_value=v, mode="discrete")
         for k, v in [("P_free", 0), ("P_rep", 1), ("M", 0), ("Y", 9), ("X", 0)]}
    m.add_species(list(s.values()))
    R = g.Reaction
    m.add_reaction([
        R(name="R1_repress", reactants={s["P_free"]: 1}, products={s["P_rep"]: 1},
          propensity_function="k_rep_max/(1+pow(X/x0,n))*P_free"),
        R(name="R2_derepress", reactants={s["P_rep"]: 1}, products={s["P_free"]: 1}, rate="k_unrep"),
        R(name="R3_transcribe", reactants={s["P_free"]: 1}, products={s["P_free"]: 1, s["M"]: 1}, rate="k_tx_eff"),
        R(name="R4_mrna_decay", reactants={s["M"]: 1}, products={}, rate="d_m"),
        R(name="R5_translate", reactants={s["M"]: 1}, products={s["M"]: 1, s["Y"]: 1}, rate="k_tl"),
        R(name="R6_lacY_dilute", reactants={s["Y"]: 1}, products={}, rate="d_y"),
        R(name="R7_uptake", reactants={s["Y"]: 1}, products={s["Y"]: 1, s["X"]: 1}, rate="k_up_eff"),
        R(name="R8_tmg_loss", reactants={s["X"]: 1}, products={}, rate="d_x"),
    ])
    # A2-validated CRP targets: mRNA made at the Hill rate set by c_crp(G), decays like lac mRNA.
    if outputs:
        tpm = crp_output_tpm(G, outputs, p)
        for gname, level in tpm.items():
            sp = g.Species(name=f"mRNA_{gname}", initial_value=0, mode="discrete")
            m.add_species(sp)
            m.add_parameter(g.Parameter(name=f"k_{gname}", expression=float(level * tpm_to_molecules * p["d_m"])))
            m.add_reaction([
                R(name=f"tx_{gname}", reactants={}, products={sp: 1}, rate=f"k_{gname}"),
                R(name=f"deg_{gname}", reactants={sp: 1}, products={}, rate="d_m"),
            ])
    return m
