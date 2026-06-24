#!/usr/bin/env python3
"""
plot_cf_cone.py

Plot the skin-friction coefficient along a blunted cone from an OpenFOAM
'surfaces' raw sample of the frictionCoefficient field.

    Cf = |tau_w,tangential| / (0.5 * rho_inf * U_inf^2)

The function object already does the normalization and normal-stripping, so
this script just reads the sampled Cf and plots it against axial distance from
the stagnation point (projection of each wall-face centre onto the body axis,
measured from the nose -- located here as the face nearest the axis).
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------- CONFIG ------------------------------------ #
# Body symmetry axis index: 0 = x, 1 = y, 2 = z
AXIS = 0

# Raw sampled frictionCoefficient file (glob; latest time by numeric value).
RAW_GLOB = "postProcessing/coneSurface/*/cf_surface_cone.raw"

# Optional digitized reference: CSV "axial_cm, Cf". None to skip.
REFERENCES = [
    ("hy2Foam_cf.csv", "-", "None", "k", "CFD: hy2Foam"),    
    ("monacoDSMC_cf.csv", "none", "^", "k", "DSMC: MONACO") 
]

OUTPUT_PNG = "cf_cone_comparison.png"
# --------------------------------------------------------------------------- #


def time_from_path(path):
    name = os.path.basename(os.path.dirname(path))
    try:
        return float(name)
    except ValueError:
        return float("nan")


def latest_raw(pattern):
    files = [f for f in glob.glob(pattern) if not np.isnan(time_from_path(f))]
    if not files:
        raise FileNotFoundError(
            f"No sampled file matched '{pattern}'. "
            "Check the surfaces function ran and the path/field/surface names."
        )
    files.sort(key=time_from_path)
    return files[-1]


def load_surface(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 4:
        raise ValueError(f"Expected >= 4 columns (x y z value) in {path}.")
    return data[:, :3], data[:, 3]


def axial_distance(coords, axis):
    """Signed axial distance from the nose (face nearest the axis), in cm."""
    x_axis = coords[:, axis]
    others = [i for i in range(3) if i != axis]
    radius = np.sqrt(coords[:, others[0]]**2 + coords[:, others[1]]**2)
    x_nose = x_axis[np.argmin(radius)]
    direction = np.sign(np.median(x_axis) - x_nose) or 1.0
    return (x_axis - x_nose) * direction * 100.0


def main():
    raw = latest_raw(RAW_GLOB)
    print(f"Reading: {raw}  (time = {time_from_path(raw):g})")

    coords, cf = load_surface(raw)
    print(f"Faces sampled: {cf.size}")
    print(f"Cf range: min = {cf.min():.4g}, max = {cf.max():.4g}")

    axial_cm = axial_distance(coords, AXIS)

    order = np.argsort(axial_cm)
    axial_cm, cf = axial_cm[order], cf[order]

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(axial_cm, cf, "-", color="r", lw=1.4, label="CFD: nexaFoam")

    for ref_file, ls, marker, color, ref_label in REFERENCES:
            if ref_file and os.path.exists(ref_file):
                ref_data = np.loadtxt(ref_file, delimiter=",")
                ax.plot(ref_data[:, 0], ref_data[:, 1], 
                        linestyle=ls, marker=marker, color=color, 
                        ms=5, lw=1.4, label=ref_label)

    ax.set_xlabel("Axial distance from stagnation point  [cm]")
    ax.set_ylabel(r"Skin-friction coefficient  $C_f$")
    ax.set_xlim(0.0, 4.0)
    ax.set_ylim(bottom=0.0)
    ax.set_title(r"non-reacting $\mathrm{N_2}$")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=150)
    print(f"Wrote {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
