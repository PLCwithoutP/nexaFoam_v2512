#!/usr/bin/env python3
"""
plot_cp_cone.py

Reproduce a Casseau-style surface pressure-coefficient distribution along a
blunted cone from an OpenFOAM 'surfaces' raw sample of the pressure field p.

    Cp = (p - p_inf) / (0.5 * rho_inf * U_inf^2)            (Casseau 2016)

"Axial distance from stagnation point" is the projection of each wall-face
centre onto the body symmetry axis, measured from the nose (stagnation) point
-- NOT arc length along the surface.

Workflow:
  1. A `surfaces` function object samples p on the cone patch (surfaceFormat raw).
  2. This script reads that raw file, builds Cp and the axial coordinate, and
     plots Cp vs axial distance, optionally overlaying digitized reference data.
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------- CONFIG ------------------------------------ #
# Freestream (Casseau 2016 Part 2, Mach 11.3 blunted cone, non-reacting N2)
P_INF   = 21.9139      # Pa
RHO_INF = 5.113e-4     # kg/m^3
U_INF   = 2764.5       # m/s

# Body symmetry axis index: 0 = x, 1 = y, 2 = z
AXIS = 0

# Path to the raw sampled pressure file (glob; the latest match is used).
# Verify the exact name written by your surfaces function object -- it is
# usually <field>_<surfaceName>.raw under postProcessing/<funcName>/<time>/.
RAW_GLOB = "postProcessing/coneSurface/*/p_cone.raw"

# Normalize the curve by its stagnation (peak) value so it peaks at 1.0, the
# way the Casseau figure is drawn. Set False to plot the absolute Cp instead
# (stagnation Cp ~ 1.84 for this case).
NORMALIZE_TO_STAGNATION = False

# Optional digitized reference: CSV with two columns "axial_cm, Cp".
# Set to None to skip. If NORMALIZE_TO_STAGNATION is True the reference is
# assumed to already be on the 0..1 scale of the published figure.
REFERENCES = [
    ("hy2Foam_cp.csv", "-", "None", "k", "CFD: hy2Foam"),    
    ("monacoDSMC_cp.csv", "none", "^", "k", "DSMC: MONACO") 
]

OUTPUT_PNG = "cp_cone_comparison.png"
# --------------------------------------------------------------------------- #


def latest_raw(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No sampled file matched '{pattern}'.\n"
            "Check that the surfaces function ran and that the path / field "
            "name / surface name in RAW_GLOB are correct."
        )
    return files[-1]


def load_surface(path, axis):
    # raw format: comment lines start with '#', data columns are  x y z value
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:          # single face -> promote to 2-D
        data = data[None, :]
    if data.shape[1] < 4:
        raise ValueError(
            f"Expected at least 4 columns (x y z p) in {path}, "
            f"got {data.shape[1]}. Did you sample a vector field by mistake?"
        )
    coords = data[:, :3]
    p = data[:, 3]
    return coords[:, axis], p


def main():
    q_inf = 0.5 * RHO_INF * U_INF**2
    print(f"q_inf = {q_inf:.6g} Pa")

    raw = latest_raw(RAW_GLOB)
    print(f"Reading: {raw}")
    x_axis, p = load_surface(raw, AXIS)
    print(f"Faces sampled: {x_axis.size}")

    # Casseau's pressure coefficient
    cp = (p - P_INF) / q_inf

    # Stagnation point = location of maximum pressure (orientation-agnostic)
    i_stag = int(np.argmax(p))
    x_stag = x_axis[i_stag]
    cp_stag = cp[i_stag]
    print(f"Stagnation: axis coord = {x_stag:.6g} m, p = {p[i_stag]:.6g} Pa, "
          f"Cp_stag = {cp_stag:.4f}")

    # Signed axial distance, positive downstream, in centimetres
    direction = np.sign(np.median(x_axis) - x_stag) or 1.0
    axial_cm = (x_axis - x_stag) * direction * 100.0

    if NORMALIZE_TO_STAGNATION:
        cp_plot = cp / cp_stag
        ylabel = r"$C_p\,/\,C_{p,\mathrm{stag}}$"
    else:
        cp_plot = cp
        ylabel = r"Pressure coefficient  $C_p$"

    order = np.argsort(axial_cm)
    axial_cm, cp_plot = axial_cm[order], cp_plot[order]

    # ------------------------------ plot ----------------------------------- #
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(axial_cm, cp_plot, "-", color="r", lw=1.4, label="CFD: nexaFoam")

    for ref_file, ls, marker, color, ref_label in REFERENCES:
            if ref_file and os.path.exists(ref_file):
                ref_data = np.loadtxt(ref_file, delimiter=",")
                ax.plot(ref_data[:, 0], ref_data[:, 1], 
                        linestyle=ls, marker=marker, color=color, 
                        ms=5, lw=1.4, label=ref_label)
                        
    ax.set_xlabel("Axial distance from stagnation point  [cm]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_title(r"non-reacting $\mathrm{N_2}$")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=150)
    print(f"Wrote {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
