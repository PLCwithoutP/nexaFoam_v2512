"""
sodAnalytical.py
----------------
Exact Riemann solver for the Sod shock tube, compared against nexaFoam output.

Run from the case root directory:
    python3 sodAnalytical.py

Input  : postProcessing/xLineSample/0.007/xLine_TTR_p_rho_U.csv
Output : postProcessing/xLineSample/0.007/sod_comparison.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, sys

# ── problem parameters ─────────────────────────────────────────────────────
GAMMA = 1.4
R_GAS = 287.0        # J/(kg·K)
T_END = 0.007        # 7 ms

RHO_L, P_L, U_L = 1.0,    1.0e5, 0.0
RHO_R, P_R, U_R = 0.125,  1.0e4, 0.0

N_EXACT = 1000       # resolution of the analytical curve

# ── paths ──────────────────────────────────────────────────────────────────
CSV_PATH = os.path.join(
    "postProcessing", "xLineSample", "0.007",
    "xLine_TTR_p_rho_U.csv"
)
OUT_PATH = os.path.join(
    "postProcessing", "xLineSample", "0.007",
    "sod_comparison.png"
)

# ══════════════════════════════════════════════════════════════════════════
#  EXACT RIEMANN SOLVER
# ══════════════════════════════════════════════════════════════════════════

def exact_riemann(x_arr, t, rho_L, p_L, u_L, rho_R, p_R, u_R, gamma, R):
    """
    Exact solution of the 1D Euler Riemann problem.
    Diaphragm at x = 0, t = 0.
    Returns rho, p, u, T arrays evaluated at positions x_arr and time t.

    Algorithm: Toro (2009), Chapter 4 — iterative pressure solver with
    Newton–Raphson, then closed-form wave sampling.
    """
    a_L = np.sqrt(gamma * p_L / rho_L)
    a_R = np.sqrt(gamma * p_R / rho_R)

    # ── pressure functions f_K and their derivatives (Toro Eqs. 4.6, 4.37) ──
    def f_K(p, p_K, rho_K, a_K):
        """Contribution of wave K to the velocity jump condition."""
        if p <= p_K:                       # rarefaction
            return (2.0*a_K/(gamma-1)) * ((p/p_K)**((gamma-1)/(2*gamma)) - 1.0)
        else:                              # shock
            A_K = 2.0 / ((gamma+1)*rho_K)
            B_K = (gamma-1)/(gamma+1) * p_K
            return (p - p_K) * np.sqrt(A_K / (p + B_K))

    def df_K(p, p_K, rho_K, a_K):
        """Derivative df_K/dp."""
        if p <= p_K:                       # rarefaction
            return (1.0/(rho_K*a_K)) * (p/p_K)**(-(gamma+1)/(2*gamma))
        else:                              # shock
            A_K = 2.0 / ((gamma+1)*rho_K)
            B_K = (gamma-1)/(gamma+1) * p_K
            sq  = np.sqrt(A_K / (p + B_K))
            return sq * (1.0 - (p - p_K) / (2.0*(p + B_K)))

    # Total pressure function whose root is p*
    f_total  = lambda p: f_K(p,p_L,rho_L,a_L) + f_K(p,p_R,rho_R,a_R) + (u_R - u_L)
    df_total = lambda p: df_K(p,p_L,rho_L,a_L) + df_K(p,p_R,rho_R,a_R)

    # ── Newton–Raphson for p* ──────────────────────────────────────────────
    # Initial guess: two-rarefaction approximation (Toro Eq. 4.46)
    p_star = ((a_L + a_R - (gamma-1)/2.0*(u_R - u_L)) /
              (a_L/p_L**((gamma-1)/(2*gamma)) +
               a_R/p_R**((gamma-1)/(2*gamma))))**(2*gamma/(gamma-1))
    p_star = max(p_star, 1e-8)

    for _ in range(100):
        fp   = f_total(p_star)
        dfp  = df_total(p_star)
        p_new = p_star - fp / dfp
        p_new = max(p_new, 1e-10)
        if abs(p_new - p_star) < 1e-10 * p_star:
            break
        p_star = p_new
    p_star = p_new

    # Star velocity
    u_star = (0.5*(u_L + u_R)
              + 0.5*(f_K(p_star,p_R,rho_R,a_R) - f_K(p_star,p_L,rho_L,a_L)))

    # ── star densities ─────────────────────────────────────────────────────
    rho_star_L = rho_L * (p_star/p_L)**(1.0/gamma)           # isentropic (rarefaction side)

    mu2 = (gamma-1)/(gamma+1)                                 # Rankine–Hugoniot (shock side)
    rho_star_R = rho_R * (p_star/p_R + mu2) / (mu2*p_star/p_R + 1.0)

    # ── wave speeds ────────────────────────────────────────────────────────
    a_star_L = a_L * (p_star/p_L)**((gamma-1)/(2*gamma))

    S_HL = u_L - a_L                            # head of left rarefaction
    S_TL = u_star - a_star_L                    # tail of left rarefaction
    S_C  = u_star                               # contact discontinuity
    S_R  = (u_R + a_R *
            np.sqrt((gamma+1)/(2*gamma) * p_star/p_R
                    + (gamma-1)/(2*gamma)))      # right-going shock

    print("── Exact Riemann solution ───────────────────────────────────")
    print(f"  p*      = {p_star:10.2f}  Pa")
    print(f"  u*      = {u_star:10.4f}  m/s")
    print(f"  rho*_L  = {rho_star_L:10.5f}  kg/m³  (behind rarefaction)")
    print(f"  rho*_R  = {rho_star_R:10.5f}  kg/m³  (behind shock)")
    print(f"  S_HL    = {S_HL:10.4f}  m/s   (rarefaction head)")
    print(f"  S_TL    = {S_TL:10.4f}  m/s   (rarefaction tail)")
    print(f"  S_C     = {S_C:10.4f}  m/s   (contact)")
    print(f"  S_R     = {S_R:10.4f}  m/s   (shock)")
    print("─────────────────────────────────────────────────────────────")

    # ── sample solution ────────────────────────────────────────────────────
    xi = x_arr / t      # self-similar variable

    rho_out = np.empty_like(x_arr)
    p_out   = np.empty_like(x_arr)
    u_out   = np.empty_like(x_arr)

    for i, s in enumerate(xi):
        if s <= S_HL:
            # Left undisturbed state
            rho_out[i], p_out[i], u_out[i] = rho_L, p_L, u_L

        elif s <= S_TL:
            # Inside rarefaction fan (Toro Eqs. 4.56–4.58)
            # Riemann invariant: u + 2a/(γ-1) = u_L + 2a_L/(γ-1)
            # Fan characteristic: ξ = u - a  →  a = u - ξ
            # Solving: u = (γ-1)/(γ+1)*u_L + 2/(γ+1)*(a_L + ξ)
            u_f   = (gamma-1)/(gamma+1)*u_L + 2.0/(gamma+1)*(a_L + s)
            a_f   = u_f - s                                # a = u - ξ
            rho_out[i] = rho_L * (a_f/a_L)**(2.0/(gamma-1))
            p_out[i]   = p_L   * (a_f/a_L)**(2.0*gamma/(gamma-1))
            u_out[i]   = u_f

        elif s <= S_C:
            # Star-left region (constant state between fan tail and contact)
            rho_out[i], p_out[i], u_out[i] = rho_star_L, p_star, u_star

        elif s <= S_R:
            # Star-right region (constant state between contact and shock)
            rho_out[i], p_out[i], u_out[i] = rho_star_R, p_star, u_star

        else:
            # Right undisturbed state
            rho_out[i], p_out[i], u_out[i] = rho_R, p_R, u_R

    T_out = p_out / (rho_out * R)
    return rho_out, p_out, u_out, T_out


# ══════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════

# ── analytical solution ────────────────────────────────────────────────────
x_exact = np.linspace(-5.0, 5.0, N_EXACT)
rho_ex, p_ex, u_ex, T_ex = exact_riemann(
    x_exact, T_END,
    RHO_L, P_L, U_L,
    RHO_R, P_R, U_R,
    GAMMA, R_GAS
)

# ── nexaFoam data ──────────────────────────────────────────────────────────
if not os.path.isfile(CSV_PATH):
    print(f"WARNING: {CSV_PATH} not found — plotting analytical solution only.")
    df = None
else:
    df = pd.read_csv(CSV_PATH)

# ── plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
fig.suptitle("Sod Shock Tube — nexaFoam vs Exact Riemann   (t = 7 ms)",
             fontsize=13)

panels = [
    (axes[0,0], p_ex,   "p",   "Pressure  $p$  [Pa]",        "tab:blue"),
    (axes[0,1], T_ex,   "TTR", "Temperature  $T_{TR}$  [K]",  "tab:orange"),
    (axes[1,0], rho_ex, "rho", "Density  $\\rho$  [kg/m³]",   "tab:green"),
    (axes[1,1], u_ex,   "U_0", "Velocity  $U_x$  [m/s]",      "tab:red"),
]

for ax, exact_data, col_name, ylabel, color in panels:
    ax.plot(x_exact, exact_data,
            color="black", linewidth=1.6, linestyle="-", label="Exact Riemann",
            zorder=3)
    if df is not None:
        ax.plot(df["x"], df[col_name],
                color=color, linewidth=1.2, linestyle="--",
                marker="o", markersize=2.5, label="nexaFoam",
                zorder=2)
    ax.axvline(0, color="gray", linewidth=0.7, linestyle=":", alpha=0.6,
               label="diaphragm x=0")
    ax.set_xlabel("$x$  [m]", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xlim(-5, 5)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)
    ax.legend(fontsize=8)

#plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150)
print(f"\nSaved → {OUT_PATH}")