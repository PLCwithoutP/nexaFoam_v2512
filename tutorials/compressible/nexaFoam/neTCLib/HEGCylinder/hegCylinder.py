#!/usr/bin/env python3
"""
cylinderSurfacePlots.py  (HEG cylinder only)

Surface pressure [kPa] and surface heat flux [MW/m^2] versus circumferential
angle theta for the nexaFoam HEG air cylinder, overlaid with SU2-NEMO and
experimental reference data (Maier 2021 / Karl).

Geometry:
  - theta is measured from the windward stagnation point, defined purely from
    the cylinder geometry and the known freestream direction (NOT from the
    pressure field), so a corner blow-up can never hijack theta=0.
  - The cylinder centre is recovered by an algebraic (Kasa) circle fit.
  - The case is now run on a HALF cylinder (theta 0..180), but only the
    windward 0..90 window is plotted, to match the reference data extent.

Styling (fixed):
  - nexaFoam     : red solid line
  - SU2-NEMO     : dark-blue dashed line
  - Experimental : black markers (no line)

Usage
    python cylinderSurfacePlots.py
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------- CONFIG ------------------------------------ #
FUNC_DIR   = "cylinderSurface"        # name of the `surfaces` function object
SURF       = "cylinder"               # surface name inside that FO
RAW_GLOB   = f"postProcessing/{FUNC_DIR}/*"

FLOW_DIR   = (1.0, 0.0, 0.0)          # freestream velocity direction (mesh frame)
THETA_MAX  = 90.0                     # plot window [deg]; data beyond is clipped
TITLE      = "HEG cylinder  (reacting air)"

# nexaFoam line style
NEXA_KW = dict(color="red", linestyle="-", lw=1.6, label="nexaFoam")

# Per-quantity definition:
#   field    : raw filename stem -> <field>_<surf>.raw
#   scale    : multiply nexaFoam values (unit conversion; refs are already in
#              the target unit and are NOT scaled)
#   take_abs : abs() the nexaFoam values (wall-heat-flux sign convention)
#   refs     : list of (csv_file, style_dict)  -- plotted if the file exists
QUANTITIES = [
    dict(
        key="pressure", field="p", scale=1e-3, take_abs=False,
        ylabel=r"Surface pressure  [kPa]",
        refs=[
            ("su2nemo_p.csv",
             dict(color="darkblue", linestyle="--", lw=1.6, label="SU2-NEMO")),
            ("experimental_p.csv",
             dict(color="k", linestyle="none", marker="o", ms=6, label="Experimental")),
        ],
    ),
    dict(
        key="heatflux", field="wallHeatFlux", scale=1e-6, take_abs=True,
        ylabel=r"Surface heat flux  [MW/m$^2$]",
        refs=[
            ("su2nemo_q.csv",
             dict(color="darkblue", linestyle="--", lw=1.6, label="SU2-NEMO")),
            ("experimental_q.csv",
             dict(color="k", linestyle="none", marker="o", ms=6, label="Experimental")),
        ],
    ),
]
# --------------------------------------------------------------------------- #


# ------------------------------ I/O helpers -------------------------------- #
def time_from_path(path):
    name = os.path.basename(os.path.dirname(path))
    try:
        return float(name)
    except ValueError:
        return float("nan")


def latest_raw(field):
    """Newest <field>_<SURF>.raw under postProcessing/<FUNC_DIR>/<time>/."""
    pattern = os.path.join(RAW_GLOB, f"{field}_{SURF}.raw")
    files = [f for f in glob.glob(pattern) if not np.isnan(time_from_path(f))]
    if not files:
        raise FileNotFoundError(
            f"No sampled file matched '{pattern}'. "
            "Check the surfaces FO ran and that the func/field/surface names match."
        )
    files.sort(key=time_from_path)
    return files[-1]


def load_surface(path):
    """raw scalar surface: columns x y z value -> (coords[N,3], value[N])."""
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 4:
        raise ValueError(
            f"Expected >= 4 columns (x y z value) in {path}, got {data.shape[1]}."
        )
    return data[:, :3], data[:, 3]


# --------------------------- geometry / theta ------------------------------ #
def in_plane_indices(coords):
    """Axis = thinnest extent (empty/wedge direction); return the other two."""
    span = coords.max(axis=0) - coords.min(axis=0)
    axis = int(np.argmin(span))
    return [i for i in range(3) if i != axis]


def circle_center(ip):
    """Algebraic (Kasa) circle fit: x^2+y^2 = a*x + b*y + c -> centre (a/2,b/2)."""
    x, y = ip[:, 0], ip[:, 1]
    A = np.column_stack((x, y, np.ones_like(x)))
    a, b, _ = np.linalg.lstsq(A, x**2 + y**2, rcond=None)[0]
    return np.array([0.5 * a, 0.5 * b])


def theta_from_upstream(coords, ip_idx, centre, flow_ip):
    """Angle (deg, 0..180) between each face's radial direction and the
    UPSTREAM direction (-flow). 0 = windward stagnation, 180 = leeward."""
    r = coords[:, ip_idx] - centre
    rn = np.linalg.norm(r, axis=1)
    rn[rn == 0.0] = 1.0
    rhat = r / rn[:, None]
    upstream = -flow_ip / np.linalg.norm(flow_ip)
    return np.degrees(np.arccos(np.clip(rhat @ upstream, -1.0, 1.0)))


# --------------------------------- driver ---------------------------------- #
def main():
    # geometry reference from the pressure sample's coordinates (never corrupt)
    pcoords, pvals = load_surface(latest_raw("p"))
    ip_idx = in_plane_indices(pcoords)
    centre = circle_center(pcoords[:, ip_idx])
    flow_ip = np.asarray(FLOW_DIR, float)[ip_idx]
    if np.linalg.norm(flow_ip) < 1e-12:
        raise SystemExit("FLOW_DIR has no in-plane component; set it correctly.")

    th_p = theta_from_upstream(pcoords, ip_idx, centre, flow_ip)
    i_pmax = int(np.argmax(pvals))
    print(f"centre = ({centre[0]:.4g}, {centre[1]:.4g}),  "
          f"theta=0 from freestream {FLOW_DIR}")
    print(f"theta span = [{th_p.min():.1f}, {th_p.max():.1f}] deg "
          f"(plotting 0..{THETA_MAX:.0f})")
    print(f"max pressure {pvals.max():.4g} Pa at theta = {th_p[i_pmax]:.1f} deg  "
          + ("OK" if th_p[i_pmax] < 10 else
             "WARNING: peak away from stagnation (corner/blow-up?)"))

    for q in QUANTITIES:
        try:
            raw = latest_raw(q["field"])
        except FileNotFoundError as e:
            print(f"[skip {q['key']}] {e}")
            continue

        coords, val = load_surface(raw)
        th = theta_from_upstream(coords, ip_idx, centre, flow_ip)
        y = np.abs(val) * q["scale"] if q["take_abs"] else val * q["scale"]

        # windward window only, sorted
        m = th <= THETA_MAX
        th_w, y_w = th[m], y[m]
        o = np.argsort(th_w)
        th_w, y_w = th_w[o], y_w[o]
        print(f"{q['key']:9s} <- {os.path.basename(raw)}  "
              f"(N={val.size}, in 0..{THETA_MAX:.0f}: {m.sum()} faces, "
              f"range [{y_w.min():.4g}, {y_w.max():.4g}])")

        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        ax.plot(th_w, y_w, **NEXA_KW)

        for ref_file, style in q["refs"]:
            if os.path.exists(ref_file):
                ref = np.loadtxt(ref_file, delimiter=",")
                if ref.ndim == 1:
                    ref = ref[None, :]
                ax.plot(ref[:, 0], ref[:, 1], **style)
            else:
                print(f"   [ref missing] {ref_file}")

        ax.set_xlabel(r"$\theta$ from stagnation point  [deg]")
        ax.set_ylabel(q["ylabel"])
        ax.set_xlim(0.0, THETA_MAX)
        ax.relim()
        ax.autoscale(axis="y", tight=False)   # y follows data (+ refs)
        ax.margins(y=0.05)
        ax.set_title(TITLE)
        ax.legend(frameon=True)
        fig.tight_layout()

        out = f"heg_{q['key']}_theta.png"
        fig.savefig(out, dpi=150)
        print(f"          -> wrote {out}")
        plt.close(fig)


if __name__ == "__main__":
    main()
