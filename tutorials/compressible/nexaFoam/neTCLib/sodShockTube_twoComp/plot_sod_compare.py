import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---- EDIT THESE ----
time = "0.0559795"
base = f"postProcessing/sampleDict/{time}"
p_file = f"{base}/xLine_p.csv"
u_file = f"{base}/xLine_U.csv"

# Normalization (optional, set to 1.0 if you want dimensional plots)
p_ref = 1.0      # e.g. left pressure pL
u_ref = 1.0      # e.g. left sound speed aL
normalize_x_to_01 = True  # make x go from 0 to 1 like the paper
# --------------------

# Load
p_df = pd.read_csv(p_file)
u_df = pd.read_csv(u_file)

# Identify columns
x_p = p_df.iloc[:, 0].to_numpy()
p   = p_df.iloc[:, 1].to_numpy()

x_u = u_df.iloc[:, 0].to_numpy()

# Find the x-velocity column automatically (first non-x column)
u_cols = [c for c in u_df.columns if c != u_df.columns[0]]
Ux = u_df[u_cols[0]].to_numpy()  # assumes first component is along x

# Sort by x (sampling sometimes already sorted, but be safe)
ip = np.argsort(x_p)
iu = np.argsort(x_u)

x_p, p  = x_p[ip], p[ip]
x_u, Ux = x_u[iu], Ux[iu]

# Optional: map x from [-5,5] to [0,1] to match the paper axis
if normalize_x_to_01:
    xmin = min(x_p.min(), x_u.min())
    xmax = max(x_p.max(), x_u.max())
    x_p = (x_p - xmin) / (xmax - xmin)
    x_u = (x_u - xmin) / (xmax - xmin)

# Normalize p and u if desired
p_plot  = p  / p_ref
u_plot  = Ux / u_ref

# Plot (paper-like: 2 panels, faint grid, legends)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.2), sharex=True)

ax[0].plot(x_p, p_plot, label="OpenFOAM")
ax[0].set_title("Sod shock test - pressure")
ax[0].set_xlabel("x")
ax[0].set_ylabel("p" if (p_ref == 1.0) else "p*")
ax[0].grid(True, alpha=0.35)

ax[1].plot(x_u, u_plot, label="OpenFOAM")
ax[1].set_title("Sod shock test - velocity")
ax[1].set_xlabel("x")
ax[1].set_ylabel("Ux" if (u_ref == 1.0) else "u*")
ax[1].grid(True, alpha=0.35)

ax[0].legend()
ax[1].legend()

plt.tight_layout()
plt.savefig("sod_two_panel.png", dpi=200)
plt.show()

