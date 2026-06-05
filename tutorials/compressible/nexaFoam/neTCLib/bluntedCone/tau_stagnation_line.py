#!/usr/bin/env python3
"""
plot_tau_cone.py

Plot the wall shear-stress magnitude along a blunted cone from an OpenFOAM
'surfaces' raw sample of the wallShearStress (vector) field.

The sampled field is the wall traction tau_w = n . Reff [Pa]; this script
plots |tau_w| against axial distance from the stagnation point (projection of
each wall-face centre onto the body axis, measured from the nose -- located
here as the face nearest the axis).
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------- CONFIG ------------------------------------ #
AXIS = 0   # body symmetry axis index: 0 = x, 1 = y, 2 = z

# Raw sampled wallShearStress file (glob; latest time by numeric value).
RAW_GLOB = "postProcessing/coneSurface/*/tau_surface_cone.raw"

REFERENCE_CSV = None       # optional "axial_cm, tau_Pa"
OUTPUT_PNG = "tau_cone.png"
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


def load_vector_surface(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 6:
        raise ValueError(
            f"Expected 6 columns (x y z vx vy vz) in {path}, "
            f"got {data.shape[1]}. Is this a vector-field sample?"
        )
    coords = data[:, :3]
    vec = data[:, 3:6]
    return coords, np.linalg.norm(vec, axis=1)


def axial_distance(coords, axis):
    x_axis = coords[:, axis]
    others = [i for i in range(3) if i != axis]
    radius = np.sqrt(coords[:, others[0]]**2 + coords[:, others[1]]**2)
    x_nose = x_axis[np.argmin(radius)]
    direction = np.sign(np.median(x_axis) - x_nose) or 1.0
    return (x_axis - x_nose) * direction * 100.0


def main():
    raw = latest_raw(RAW_GLOB)
    print(f"Reading: {raw}  (time = {time_from_path(raw):g})")

    coords, tau = load_vector_surface(raw)
    print(f"Faces sampled: {tau.size}")
    print(f"|tau_w| range: min = {tau.min():.4g} Pa, max = {tau.max():.4g} Pa")

    axial_cm = axial_distance(coords, AXIS)

    order = np.argsort(axial_cm)
    axial_cm, tau = axial_cm[order], tau[order]

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(axial_cm, tau, "-", color="k", lw=1.4, label="nexaFoam")

    if REFERENCE_CSV and os.path.exists(REFERENCE_CSV):
        ref = np.loadtxt(REFERENCE_CSV, delimiter=",")
        ax.plot(ref[:, 0], ref[:, 1], "^", color="k", ms=5, ls="none",
                label="reference")

    ax.set_xlabel("Axial distance from stagnation point  [cm]")
    ax.set_ylabel(r"Wall shear stress  $|\tau_w|$  [Pa]")
    ax.set_xlim(0.0, float(axial_cm.max()))
    ax.set_ylim(bottom=0.0)
    ax.set_title(r"non-reacting $\mathrm{N_2}$")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=150)
    print(f"Wrote {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
