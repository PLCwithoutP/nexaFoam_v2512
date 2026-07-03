#!/usr/bin/env python3
"""
Cone-surface wall heat-flux distribution vs surface distance, several time
snapshots, in the style of the Vatansever (2020) / hyperReactingFoam double-cone
validation.

Reads OpenFOAM sampled-surface 'raw' output for the cone patch (wallHeatFlux)
and overlays time snapshots so you can watch convergence toward steady state.

Usage:
    python q_time_plot.py
Edit the CONFIG block for your paths / surface name / times.
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
TIMES       = None      # None = all available; or e.g. [1.17e-4, 2.83e-4, 1.75e-3]
TIME_RTOL   = 0.05      # relative tolerance when matching requested times
OUT_PNG     = "time_dependent_q.png"
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


def pick_times(all_times):
    if TIMES is None:
        return all_times
    chosen = []
    for treq in TIMES:
        cand = min(all_times, key=lambda td: abs(td[0] - treq))
        if abs(cand[0] - treq) <= TIME_RTOL * max(abs(treq), 1e-30):
            chosen.append(cand)
        else:
            print(f"  [warn] no snapshot near t={treq:g}s (closest {cand[0]:g})")
    return chosen


def main():
    all_times = available_times()
    if not all_times:
        raise SystemExit(f"No time directories under {PP_DIR!r}")
    times = pick_times(all_times)

    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.get_cmap("turbo")
    n = max(len(times), 1)

    for i, (t, tdir) in enumerate(times):
        path = find_raw_file(tdir)
        if path is None:
            print(f"  [skip] no {FIELD} file in {tdir}")
            continue
        x, y, _, q = load_surface(path)
        order, s = surface_coordinate(x, y)
        ax.plot(s * LEN_TO_CM + S_OFFSET_CM, Q_SIGN * q[order] * Q_TO_WCM2,
                color=cmap(i / n), lw=1.6,
                label=f"nexaFoam : qw ~ {t:.2e}s")

    ax.set_xlabel("Horizontal distance on cone surface (cm)", fontweight="bold")
    ax.set_ylabel("Heat Flux (W/cm$^2$)", fontweight="bold")
    ax.set_xlim(4.0, 16.0)
    ax.set_ylim(0,1500)
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=9, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
