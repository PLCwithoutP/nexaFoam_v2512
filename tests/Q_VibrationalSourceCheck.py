import math

# Import thermodynamic / transport helper functions and constants
import thermodynamicPropertiesCheck as thermo

# ------------------------------------------------------------
# Mixture setup and global constants
# ------------------------------------------------------------

# Example N2/O2 mixture (mass fractions)
Y = {"N2": 0.45, "O2": 0.55}

# Species list
species = ["N2", "O2"]

# Molar masses [kg/mol] – from thermodynamicPropertiesCheck.py
M = {
    "N2": thermo.molWeight_N2,
    "O2": thermo.molWeight_O2,
}

# Species-specific gas constants [J/(kg K)] – already per unit mass
R_s = {
    "N2": thermo.R_N2,
    "O2": thermo.R_O2,
}

# Vibrational characteristic temperatures [K]
theta_v = {
    "N2": thermo.thetaVib_N2,
    "O2": thermo.thetaVib_O2,
}

# Universal constants
R_univ = 8.314462618  # J/(mol K)
kB     = 1.380649e-23 # J/K
pi     = math.pi

# --- Knab VV constants (matching your C++ ctor) ---
NA          = 6.022e23   # Avogadro
sigma_vv    = 3.0e-21    # σ'_m,l  (N2, O2, NO)
P_ml_const  = 1.0e-2     # P_{m,l}

# State inputs (can be randomized if desired)
T_tr = 2000.0   # translational-rotational temperature [K]
T_v  = 1000.0   # vibrational temperature [K]
p    = 1.0      # pressure [Pa]

# ------------------------------------------------------------
# Composition and density
# ------------------------------------------------------------

# Mass -> mole fractions
den_mol = sum(Y[s] / M[s] for s in species)
X = {s: (Y[s] / M[s]) / den_mol for s in species}

# Mixture density using mixture gas constant
R_mix = sum(Y[s] * R_s[s] for s in species)  # since R_s is already J/(kg K)
rho   = p / (R_mix * T_tr)
rho_s = {s: Y[s] * rho for s in species}

# Number densities (using Boltzmann)
n_tot = p / (kB * T_tr)              # total number density [1/m^3]
n_s   = {s: X[s] * n_tot for s in species}


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def reduced_mass(s, r):
    """Reduced mass μ = Ms Mr / (Ms + Mr) [kg/mol]."""
    Ms, Mr = M[s], M[r]
    return Ms * Mr / (Ms + Mr)


# ---------- (1) Millikan–White τ^{MW}_{s-r,V-T}
def A_sr(s, r):
    mu = reduced_mass(s, r)
    return 1.16 * 10.0 ** (-4.5) * math.sqrt(mu) * theta_v[s] ** (4.0 / 3.0)


def B_sr(s, r):
    mu = reduced_mass(s, r)
    return 0.015 * 10.0 ** (-0.75) * mu ** 0.25


def tau_MW(s, r, p_local, T_tr_local):
    """Millikan–White vibrational–translational relaxation time τ^MW [s]."""
    Asr = A_sr(s, r)
    Bsr = B_sr(s, r)
    exponent = Asr * (T_tr_local ** (-1.0 / 3.0) - Bsr) - 18.42
    return (1.0 / p_local) * math.exp(exponent)


# ---------- (2) Park τ^{P}_{s-r,V-T}
def sigma_v(s, T_tr_local):
    """Temperature-dependent vibrational collision cross-section σ_v,s [m^2]."""
    sigma0 = {
        "N2": 3.0e-21,
        "O2": 3.0e-21,
    }[s]
    return sigma0 * (50000.0 / T_tr_local) ** 2


def c_bar_s(s, T_tr_local):
    """Mean thermal speed c̄_s [m/s]."""
    Ms = M[s]
    return math.sqrt(8.0 * R_univ * T_tr_local / (pi * Ms))


def v_mp_s(s, T_tr_local):
    """Most probable speed v_mp [m/s] = sqrt(2 R T / M)."""
    Ms = M[s]
    return math.sqrt(2.0 * R_univ * T_tr_local / Ms)


def n_r(r, p_local, T_tr_local):
    """Number density of species r [1/m^3] using Boltzmann constant."""
    n_tot_local = p_local / (kB * T_tr_local)
    return X[r] * n_tot_local


def tau_P(s, r, p_local, T_tr_local):
    """Park correction τ^P [s]."""
    cbar = c_bar_s(s, T_tr_local)
    sig  = sigma_v(s, T_tr_local)
    n    = n_r(r, p_local, T_tr_local)
    return 1.0 / (cbar * sig * n)


# ---------- (3) Combined τ_{s-r,V-T}
def tau_sr_total(s, r, p_local, T_tr_local):
    return tau_MW(s, r, p_local, T_tr_local) + tau_P(s, r, p_local, T_tr_local)


# ---------- (4) Mode relaxation time τ_{s,V-T}
def tau_mode_s(s, p_local, T_tr_local):
    """Mode relaxation time τ_s,V-T from species-specific τ_s-r,V-T."""
    num   = sum(X[r] for r in species)
    denom = 0.0
    for r in species:
        tau_sr = tau_sr_total(s, r, p_local, T_tr_local)
        denom += X[r] / tau_sr
    return num / denom  # equivalent to 1 / Σ_r (X_r / τ_s-r)


# ---------- Vibrational thermodynamics via thermodynamicPropertiesCheck.py
def Cv_vib_mass(s, T):
    """Vibrational Cv per unit mass using imported CvVib()."""
    theta  = theta_v[s]
    R_spec = R_s[s]
    return thermo.CvVib(T, theta, R_spec)


def e_vib_mass(s, T):
    """Vibrational energy per unit mass using θ_v and R_s (harmonic oscillator)."""
    theta  = theta_v[s]
    R_spec = R_s[s]
    z      = math.exp(theta / T)
    return R_spec * theta / (z - 1.0)


# ---------- (5) Q_s,V-T (classical & modified Landau–Teller)
def Q_classical(s, p_local, T_tr_local, T_v_local):
    tau_s = tau_mode_s(s, p_local, T_tr_local)
    ev_tr = e_vib_mass(s, T_tr_local)
    ev_v  = e_vib_mass(s, T_v_local)
    return rho_s[s] * (ev_tr - ev_v) / tau_s


def Q_modified(s, p_local, T_tr_local, T_v_local):
    tau_s = tau_mode_s(s, p_local, T_tr_local)
    Cv_v  = Cv_vib_mass(s, T_tr_local)
    return (T_tr_local / T_v_local) * (T_tr_local - T_v_local) * Cv_v * rho_s[s] / tau_s


# ------------------------------------------------------------
# New: Knab VV helper functions
# ------------------------------------------------------------

def p_s(p_local, s):
    """Partial pressure of species s: p_s = X_s * p."""
    return X[s] * p_local


def rho_from_partial(p_s_local, s, T_tr_local):
    """Species density from partial pressure: ρ_s = p_s / (R_s T)."""
    return p_s_local / (R_s[s] * T_tr_local)


def Q_VV_s(m, p_local, T_tr_local, T_v_local, verbose=False):
    """
    Knab VV source term Q_{m,V-V} for species m.
    Mirrors the C++ implementation using:
      Q_m = Σ_{l≠m} N_A σ'_m,l P_m,l c̄_m sqrt(ρ_l/M_l) ρ_m
            [ e_v,m(T_tr) e_v,l(T_v) / e_v,l(T_tr) - e_v,m(T_v) ]
    """
    # partial pressure and density of active species m
    p_m   = p_s(p_local, m)
    rho_m = rho_from_partial(p_m, m, T_tr_local)

    ev_m_Ttr = e_vib_mass(m, T_tr_local)
    ev_m_Tv  = e_vib_mass(m, T_v_local)

    cbar_m = c_bar_s(m, T_tr_local)

    Qm = 0.0

    for l in species:
        if l == m:
            continue

        # partial pressure and density of collider l
        p_l   = p_s(p_local, l)
        rho_l = rho_from_partial(p_l, l, T_tr_local)

        ev_l_Ttr = e_vib_mass(l, T_tr_local)
        ev_l_Tv  = e_vib_mass(l, T_v_local)

        if abs(ev_l_Ttr) < 1e-50:
            continue

        vel_factor = cbar_m * math.sqrt(rho_l / M[l])

        energy_bracket = ev_m_Ttr * (ev_l_Tv / ev_l_Ttr) - ev_m_Tv

        contrib = (
            NA
            * sigma_vv
            * P_ml_const
            * vel_factor
            * rho_m
            * energy_bracket
        )

        Qm += contrib

        if verbose:
            print(
                f"    VV contrib m={m}, l={l}: "
                f"vel_factor={vel_factor:.3e}, "
                f"bracket={energy_bracket:.3e}, "
                f"Q_ml={contrib:.3e} W/m^3"
            )

    return Qm


# ------------------------------------------------------------
# Main: VT + VV source evaluation
# ------------------------------------------------------------
if __name__ == "__main__":
    print("============================================================")
    print("VT & VV source term sanity script (using thermodynamicPropertiesCheck.py)")
    print("------------------------------------------------------------")
    print(f"T_tr = {T_tr:.1f} K,  T_v = {T_v:.1f} K,  p = {p:.3e} Pa")
    print("Mass fractions Y:", Y)
    print("Mole fractions X:", {s: f"{X[s]:.4f}" for s in species})
    print(f"Total density rho = {rho:.3e} kg/m^3")
    print(f"Total number density n_tot = {n_tot:.3e} 1/m^3")
    print("Number densities n_s [1/m^3]:")
    for s in species:
        print(f"  {s}: {n_s[s]:.3e}")
    print("------------------------------------------------------------")

    print("Most probable speeds v_mp [m/s] at T_tr:")
    for s in species:
        vmp = v_mp_s(s, T_tr)
        print(f"  {s}: {vmp:.3e}")
    print("------------------------------------------------------------")

    for s in species:
        print(f"\n### Species s = {s} ###")
        for r in species:
            Asr = A_sr(s, r)
            Bsr = B_sr(s, r)
            tau_mw = tau_MW(s, r, p, T_tr)
            tau_p  = tau_P(s, r, p, T_tr)
            tau_sr = tau_sr_total(s, r, p, T_tr)
            print(f"  Collision with r = {r}:")
            print(f"    A_sr       = {Asr:.3e}")
            print(f"    B_sr       = {Bsr:.3e}")
            print(f"    tau_MW     = {tau_mw:.3e} s")
            print(f"    tau_P      = {tau_p:.3e} s")
            print(f"    tau_s-r    = {tau_sr:.3e} s")

        tau_s = tau_mode_s(s, p, T_tr)
        Qc    = Q_classical(s, p, T_tr, T_v)
        Qm    = Q_modified(s, p, T_tr, T_v)
        Qvv   = Q_VV_s(s, p, T_tr, T_v, verbose=True)

        print(f"\n  Mode relaxation time tau_s,V-T = {tau_s:.3e} s")
        print(f"  Q_s,V-T (classical LT)        = {Qc:.3e} W/m^3")
        print(f"  Q_s,V-T (modified LT)         = {Qm:.3e} W/m^3")
        print(f"  Q_s,V-V (Knab)                = {Qvv:.3e} W/m^3")
        print("------------------------------------------------------------")

