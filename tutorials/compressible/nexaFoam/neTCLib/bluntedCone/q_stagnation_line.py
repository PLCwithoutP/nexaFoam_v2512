#!/usr/bin/env python3
"""
plot_q_cone.py

Plot the conductive wall heat flux along a blunted cone from an OpenFOAM
'surfaces' raw sample of the wallHeatFlux field [W/m^2], against axial distance
from the stagnation point (nose located as the face nearest the axis).
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------- CONFIG ------------------------------------ #
AXIS = 0   # body symmetry axis index: 0 = x, 1 = y, 2 = z

RAW_GLOB = "postProcessing/coneSurface/*/wallHeatFlux_cone.raw"

# Plot in kW/m^2 (set 1.0 for W/m^2)
SCALE = 1.0e-3
YLABEL = r"Wall heat flux  $q_w$  [kW/m$^2$]"

REFERENCES = [
    #("hy2Foam_ref.csv", "^", "hy2Foam"),
    #("casseau_n2_ref.csv", "s", "Casseau 2016")
]

OUTPUT_PNG = "q_cone.png"
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
    x_axis = coords[:, axis]
    others = [i for i in range(3) if i != axis]
    radius = np.sqrt(coords[:, others[0]]**2 + coords[:, others[1]]**2)
    x_nose = x_axis[np.argmin(radius)]
    direction = np.sign(np.median(x_axis) - x_nose) or 1.0
    return (x_axis - x_nose) * direction * 100.0


def main():
    raw = latest_raw(RAW_GLOB)
    print(f"Reading: {raw}  (time = {time_from_path(raw):g})")

    coords, q = load_surface(raw)
    print(f"Faces sampled: {q.size}")
    print(f"q_w range: min = {q.min():.4g} W/m^2, max = {q.max():.4g} W/m^2")

    axial_cm = axial_distance(coords, AXIS)

    order = np.argsort(axial_cm)
    axial_cm, q = axial_cm[order], q[order]*SCALE

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(axial_cm, q, "-", color="k", lw=1.4, label="nexaFoam")

    for ref_file, marker, ref_label in REFERENCES:
            if ref_file and os.path.exists(ref_file):
                ref_data = np.loadtxt(ref_file, delimiter=",")
                ax.plot(ref_data[:, 0], ref_data[:, 1], marker, color="k", ms=5, ls="none",
                        label=ref_label)

    ax.set_xlabel("Axial distance from stagnation point  [cm]")
    ax.set_ylabel(YLABEL)
    ax.set_xlim(0.0, float(axial_cm.max()))
    ax.set_title(r"non-reacting $\mathrm{N_2}$")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=150)
    print(f"Wrote {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
