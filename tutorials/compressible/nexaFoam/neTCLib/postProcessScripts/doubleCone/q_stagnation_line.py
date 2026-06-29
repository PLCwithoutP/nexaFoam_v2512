#!/usr/bin/env python3
"""
Cone-surface wall heat-flux distribution at the latest written time, compared
against digitized reference curves (experiment + reference CFD), in the style of
the Vatansever (2020) / hyperReactingFoam double-cone validation.

Reads OpenFOAM sampled-surface 'raw' output for the cone patch (wallHeatFlux).

Usage:
    python q_stagnation_line.py
Edit the CONFIG block for your paths / surface name / reference files.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------- CONFIG -------------------------------------
PP_DIR      = "postProcessing/coneSurface"  # function-object output directory
SURFACE     = "cone"                        # sampled surface name
FIELD       = "wallHeatFlux"                # field name in the raw file
Q_TO_WCM2   = 1.0e-4                         # W/m^2 -> W/cm^2  (1 m^2 = 1e4 cm^2)
Q_SIGN      = 1.0                            # set -1.0 if your q_w comes out negative
LEN_TO_CM   = 1.0e2                          # m  -> cm
X_AXIS      = "axial"                        # "arclength" (s from tip) or "axial" (x)
S_OFFSET_CM = 0.0                            # shift to align with the reference origin
# Reference curves. x assumed in CENTIMETRES, y in W/cm^2.
REF_CSVS = [
    ("hyperReactingFoam_q.csv", "Vatansever (2020), CFD"),  # 0: blue line + markers
    ("holden_q.csv",             "Holden (2014), Experiment"),  # 1: black square markers
]
REF_X_OFFSET_M = 0.004527                    # +4.527 mm added to every reference x
OUT_PNG     = "validation_graph_q.png"
# ----------------------------------------------------------------------------


def find_raw_file(time_dir):
    """Locate the raw surface file for FIELD on SURFACE in a time directory."""
    candidates = [
        os.path.join(time_dir, f"{FIELD}_{SURFACE}.raw"),
        os.path.join(time_dir, f"{SURFACE}_{FIELD}.raw"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    hits = glob.glob(os.path.join(time_dir, f"*{FIELD}*.raw"))
    return hits[0] if hits else None


def load_surface(path):
    """Return x, y, z, q from an OpenFOAM raw scalar-surface file."""
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def surface_coordinate(x, y):
    """Sort streamwise (tip = smallest x); return sort order and 1-D coord [m]."""
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    if X_AXIS == "axial":
        s = xs
    else:  # arc length from the tip
        ds = np.hypot(np.diff(xs), np.diff(ys))
        s = np.concatenate(([0.0], np.cumsum(ds)))
    return order, s


def available_times():
    out = []
    for d in glob.glob(os.path.join(PP_DIR, "*")):
        try:
            out.append((float(os.path.basename(d)), d))
        except ValueError:
            pass
    return sorted(out, key=lambda t: t[0])


def main():
    all_times = available_times()
    if not all_times:
        raise SystemExit(f"No time directories under {PP_DIR!r}")
    t, tdir = all_times[-1]            # latest written time only

    fig, ax = plt.subplots(figsize=(9, 6))

    path = find_raw_file(tdir)
    if path is None:
        raise SystemExit(f"No {FIELD} file in {tdir}")
    x, y, _, q = load_surface(path)
    order, s = surface_coordinate(x, y)
    ax.plot(s * LEN_TO_CM + S_OFFSET_CM, Q_SIGN * q[order] * Q_TO_WCM2,
            color="red", lw=2.0, label=f"nexaFoam : $q_w$ ~ {t:.2e}s")

    # per-curve styles (order matches REF_CSVS)
    ref_styles = [
        dict(color="blue", ls="-", lw=1.0, marker="+", markersize=6,
             markeredgewidth=1.2, markevery=2),
        dict(color="black", ls="none", marker="s", markersize=7,
             markerfacecolor="black"),
    ]
    for j, (rpath, rlabel) in enumerate(REF_CSVS):
        if not os.path.isfile(rpath):
            print(f"  [warn] reference not found: {rpath}")
            continue
        ref = np.loadtxt(rpath, delimiter=",")
        if ref.ndim == 1:
            ref = ref[None, :]
        xref = ref[:, 0] + REF_X_OFFSET_M * LEN_TO_CM
        ax.plot(xref, ref[:, 1], label=rlabel, **ref_styles[j % len(ref_styles)])

    ax.set_xlabel("Horizontal distance on cone surface (cm)", fontweight="bold")
    ax.set_ylabel("Heat Flux (W/cm$^2$)", fontweight="bold")
    ax.set_xlim(4.0, 16.0)
    ax.set_ylim(0,600)
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=9, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
