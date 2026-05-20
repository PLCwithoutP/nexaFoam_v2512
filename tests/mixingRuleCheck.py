import math

# Given mass fractions
Y = {"N2": 0.45, "O2": 0.55}

# Molar masses in kg/mol
M = {"N2": 28.0134e-3, "O2": 31.9988e-3}

# Dynamic viscosities of pure species at ~1000 K (Pa·s)
mu = {"N2": 4.01105e-05, "O2": 4.62162e-05}

# Thermal conductivities of pure species at ~1000 K (W/m/K)
k = {"N2": 0.0991814, "O2": 0.107939}

species = ["N2", "O2"]

# --- 1) Convert mass fractions to mole fractions ---
den = sum(Y[s] / M[s] for s in species)
x = {s: (Y[s] / M[s]) / den for s in species}

# --- 2) Wilke phi_ij function (uses viscosities and molar masses) ---
def phi(i, j):
    if i == j:
        return 1.0
    mu_i, mu_j = mu[i], mu[j]
    M_i, M_j = M[i], M[j]
    return (1.0 + math.sqrt(mu_i / mu_j) * (M_j / M_i) ** 0.25) ** 2 / \
           (math.sqrt(8.0) * math.sqrt(1.0 + M_i / M_j))

# --- 3) Wilke mixture viscosity μ_mix ---
mu_mix = 0.0
for i in species:
    denom_i = sum(x[j] * phi(i, j) for j in species)
    mu_mix += x[i] * mu[i] / denom_i

# --- 4) Wilke mixture thermal conductivity k_mix ---
k_mix = 0.0
for i in species:
    denom_i = sum(x[j] * phi(i, j) for j in species)
    k_mix += x[i] * k[i] / denom_i

print("Mole fractions:", x)
print("**************************************************************")
print("Specie viscosities:", mu)
print(f"Mixture viscosity (Wilke)          = {mu_mix:.3e} Pa·s")
print("**************************************************************")
print("Specie thermal conductivities:", k)
print(f"Mixture thermal conductivity (Wilke) = {k_mix:.5f} W/m/K")

