import math
import thermodynamicPropertiesCheck as thermo

# ------------------------------------------------------------
# Mixture setup and global constants
# ------------------------------------------------------------

# Example N2/O2 mixture (mass fractions)
Y = {"N2": 0.45, "O2": 0.55}
species = ["N2", "O2"]

# Molar masses [kg/mol] – from thermodynamicPropertiesCheck.py
M = {
    "N2": thermo.molWeight_N2,
    "O2": thermo.molWeight_O2,
}

# Species-specific gas constants [J/(kg K)] – per unit mass
R_s = {
    "N2": thermo.R_N2,
    "O2": thermo.R_O2,
}

# State
T_tr = 2000.0         # translational–rotational temperature [K]
p    = 1              # mixture pressure [Pa]

# Species Lewis numbers (you can tune these)
Le_s = {
    "N2": 1.4,
    "O2": 1.4,
}


# ------------------------------------------------------------
# Composition helpers
# ------------------------------------------------------------

def mass_to_mole_fractions(Y_local, M_local):
    """Convert mass fractions to mole fractions."""
    den = sum(Y_local[s] / M_local[s] for s in species)
    return {s: (Y_local[s] / M_local[s]) / den for s in species}


X = mass_to_mole_fractions(Y, M)


# ------------------------------------------------------------
# Thermodynamic helpers (Cp, kappa_TR, mu)
# ------------------------------------------------------------

def Cp_TR(specie, T):
    """
    Translational-rotational Cp per unit mass [J/(kg K)] for a given species.
    Uses JANAF Cp/R from thermodynamicPropertiesCheck.py.
    """
    coeffs = thermo.JanafCoeffs(T, specie)
    Cp_over_R = thermo.JanafCp(coeffs, T)   # Cp/R (dimensionless)
    R_spec = R_s[specie]
    return Cp_over_R * R_spec


def kappa_TR(specie, T):
    """
    Translational-rotational thermal conductivity k_tr [W/(m K)]
    for a given species (Blottner-Eucken-like). Reuses kappaTR()
    from thermodynamicPropertiesCheck.py.
    """
    R_spec = R_s[specie]
    Cp = Cp_TR(specie, T)
    return thermo.kappaTR(T, specie, R_spec, Cp)


def mu_specie(specie, T):
    """Dynamic viscosity μ(T) via Blottner fit from thermodynamicPropertiesCheck.py."""
    return thermo.mu(T, specie)


# ------------------------------------------------------------
# Wilke mixing rule for viscosity / conductivity
# (same structure as mixingRuleCheck.py)
# ------------------------------------------------------------

def phi(i, j, T):
    """Wilke phi_ij(T) using viscosities and molar masses."""
    if i == j:
        return 1.0
    mu_i = mu_specie(i, T)
    mu_j = mu_specie(j, T)
    M_i = M[i]
    M_j = M[j]
    return (1.0 + math.sqrt(mu_i / mu_j) * (M_j / M_i) ** 0.25) ** 2 / \
           (math.sqrt(8.0) * math.sqrt(1.0 + M_i / M_j))


def wilke_mixture_value(pure_prop_dict, T):
    """
    Generic Wilke mixer.
    pure_prop_dict: dict(specie -> pure-property value)
    Returns mixture property (for kappa_tr here).
    """
    mix_val = 0.0
    for i in species:
        denom_i = sum(X[j] * phi(i, j, T) for j in species)
        mix_val += X[i] * pure_prop_dict[i] / denom_i
    return mix_val


def mixture_kappa_TR(T):
    k_pure = {s: kappa_TR(s, T) for s in species}
    return wilke_mixture_value(k_pure, T)


# ------------------------------------------------------------
# Mixture Cp, R_mix, rho
# ------------------------------------------------------------

def mixture_Cp(T):
    """Mass-weighted mixture Cp [J/(kg K)]."""
    Cp_s = {s: Cp_TR(s, T) for s in species}
    Cp_mix = sum(Y[s] * Cp_s[s] for s in species)
    return Cp_mix, Cp_s


def mixture_R():
    """Mass-weighted mixture gas constant R_mix [J/(kg K)]."""
    return sum(Y[s] * R_s[s] for s in species)


def mixture_rho(p_local, T):
    """Ideal-gas mixture density ρ = p / (R_mix T)."""
    R_mix = mixture_R()
    return p_local / (R_mix * T)


# ------------------------------------------------------------
# Effective diffusion coefficient
# ------------------------------------------------------------

def effective_diffusion_coeff(k_tr, rho, Cp, Le):
    """
    D_s = k_tr * Le / (rho * Cp)
    Returns D_s [m^2/s] for a given Lewis number Le.
    """
    return k_tr * Le / (rho * Cp)


# ------------------------------------------------------------
# Main sanity check
# ------------------------------------------------------------

if __name__ == "__main__":
    print("============================================================")
    print("Mass diffusion sanity script: effective diffusion coefficient")
    print("D_s = D = k_tr * Le / (rho * Cp)")
    print("------------------------------------------------------------")
    print(f"T_tr = {T_tr:.1f} K,  p = {p:.3e} Pa")
    print("Mass fractions Y:", Y)
    print("Mole fractions X:", {s: f"{X[s]:.4f}" for s in species})

    # Species Cp and kappa_tr
    Cp_mix, Cp_s = mixture_Cp(T_tr)
    print("------------------------------------------------------------")
    print("Species Cp_TR and kappa_TR at T_tr:")
    for s in species:
        cp_s = Cp_s[s]
        k_s = kappa_TR(s, T_tr)
        mu_s = mu_specie(s, T_tr)
        print(f"  {s}: Cp_TR = {cp_s:10.3f} J/(kg K), "
              f"k_tr = {k_s:10.3e} W/(m K), "
              f"mu = {mu_s:10.3e} Pa·s")

    # Mixture properties
    k_tr_mix = mixture_kappa_TR(T_tr)
    rho_mix  = mixture_rho(p, T_tr)
    R_mix    = mixture_R()

    print("------------------------------------------------------------")
    print(f"Mixture Cp_TR   = {Cp_mix:.3f} J/(kg K)")
    print(f"Mixture k_tr    = {k_tr_mix:.3e} W/(m K)")
    print(f"Mixture R_mix   = {R_mix:.3f} J/(kg K)")
    print(f"Mixture density = {rho_mix:.3e} kg/m^3")

    # Effective diffusion coefficients per species (same k_tr, species Le_s)
    print("------------------------------------------------------------")
    print("Effective diffusion coefficients D_s [m^2/s]:")
    for s in species:
        Le = Le_s[s]
        D_s = effective_diffusion_coeff(k_tr_mix, rho_mix, Cp_mix, Le)
        print(f"  {s}: Le = {Le:5.2f}  ->  D_s = {D_s:.3e} m^2/s")

    print("============================================================")
    print("Check finished!")

