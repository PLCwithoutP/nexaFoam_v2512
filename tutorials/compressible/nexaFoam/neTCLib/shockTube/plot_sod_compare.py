"""
plotSodTube.py
--------------
Plots p, TTR, rho, U_0 from nexaFoam Sod shock tube results.
Run from the case root directory:
    python3 plotSodTube.py
Output: postProcessing/xLineSample/0.007/sodTube_results.png
"""

import matplotlib
matplotlib.use("Agg")          # headless — no X11 needed on WSL
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

# ── paths ──────────────────────────────────────────────────────────────────
CSV_PATH = os.path.join(
    "postProcessing", "xLineSample", "0.007",
    "xLine_TTR_p_rho_U.csv"
)
OUT_PATH = os.path.join(
    "postProcessing", "xLineSample", "0.007",
    "sodTube_results.png"
)

if not os.path.isfile(CSV_PATH):
    sys.exit(f"ERROR: cannot find {CSV_PATH}\n"
             "Run 'nexaFoam -postProcess -time 0.007' first.")

# ── load ───────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

x    = df["x"]
p    = df["p"]
TTR  = df["TTR"]
rho  = df["rho"]
U0   = df["U_0"]

# ── plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
fig.suptitle("Sod Shock Tube  —  nexaFoam   (t = 7 ms)", fontsize=13)

panels = [
    (axes[0, 0], p,   "Pressure  $p$  [Pa]",          "tab:blue"),
    (axes[0, 1], TTR, "Temperature  $T_{TR}$  [K]",   "tab:orange"),
    (axes[1, 0], rho, "Density  $\\rho$  [kg/m³]",    "tab:green"),
    (axes[1, 1], U0,  "Velocity  $U_x$  [m/s]",       "tab:red"),
]

for ax, data, ylabel, color in panels:
    ax.plot(x, data, color=color, linewidth=1.4, label="nexaFoam")
    ax.set_xlabel("$x$  [m]", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xlim(-5, 5)
    ax.axvline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.4,
               label="diaphragm")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=8)

#plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150)
print(f"Saved → {OUT_PATH}")