#!/usr/bin/env python3
"""
cylinder_time_surface_plots.py

Time-dependent surface distributions for the nexaFoam Mach-20 N2 half-cylinder,
in the convergence-watching style of the double-cone p/q snapshots: several time
snapshots overlaid with a turbo colourmap so you can watch the surface march
toward steady state.

Three figures are produced, each vs. circumferential angle theta:
    * pressure coefficient   Cp
    * skin-friction coeff.   Cf
    * surface heat flux      [W/cm^2]

theta runs 0 deg (windward stagnation) -> 180 deg (leeward). The cylinder centre
is recovered by an algebraic (Kasa) circle fit and theta is built purely from
geometry + the known freestream direction, so it is immune to spurious peaks in
the solution (e.g. a corner blow-up that would hijack a max-pressure anchor).

Usage:
    python cylinder_time_surface_plots.py
Edit the CONFIG block for your paths / surface name / times.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# =============================== CONFIG ==================================== #
RHO_INF   = 1.363e-5        # kg/m^3   (Table 2)
U_INF     = 6047.0          # m/s
FLOW_DIR  = (1.0, 0.0, 0.0)  # freestream velocity direction in the mesh frame
THETA_MAX = 180.0           # half cylinder

FUNC_DIR  = "coneSurface"  # name of the `surfaces` function object
                               # -> postProcessing/<FUNC_DIR>/<time>/<field>_<surf>.raw
SURF      = "cylinder"         # surface name inside the `surfaces` FO
PP_DIR    = f"postProcessing/{FUNC_DIR}"

TIMES     = None            # None = all available; or e.g. [1.5e-4, 2.7e-4, 4.9e-4]
TIME_RTOL = 0.05            # relative tolerance when matching requested times

Q_SIGN    = 1.0             # set -1.0 if wallHeatFlux comes out negative into the wall

# Each quantity: key (png name), field (raw stem), tag (legend), scale, take_abs,
# ylabel, and an optional fixed ylim (None -> autoscale to the data).
QUANTITIES = [
    dict(key="cp", field="cp_surface", tag="Cp", scale=1.0, take_abs=False,
         ylabel=r"Pressure coefficient  $C_p$", ylim=None),
    dict(key="cf", field="cf_surface", tag="Cf", scale=1.0, take_abs=False,
         ylabel=r"Skin-friction coefficient  $C_f$", ylim=None),
    dict(key="heatflux", field="wallHeatFlux", tag="qw", scale=1e-4, take_abs=True,
         ylabel=r"Surface heat flux  [W/cm$^2$]", ylim=None),
]
# =========================================================================== #


# ------------------------------ I/O helpers -------------------------------- #
def available_times():
    """Sorted list of (time, time_dir) under postProcessing/<func>/."""
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


def find_raw_file(time_dir, field):
    """Locate the raw surface file for `field` on SURF in a time directory."""
    candidates = [
        os.path.join(time_dir, f"{field}_{SURF}.raw"),
        os.path.join(time_dir, f"{SURF}_{field}.raw"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    hits = glob.glob(os.path.join(time_dir, f"*{field}*.raw"))
    return hits[0] if hits else None


def load_surface(path):
    """raw scalar surface: columns x y z value -> (coords[N,3], value[N])."""
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 4:
        raise ValueError(
            f"Expected >= 4 columns (x y z value) in {path}, got {data.shape[1]}. "
            "Did the FO write a vector field here?"
        )
    return data[:, :3], data[:, 3]


# --------------------------- geometry / theta ------------------------------ #
def in_plane_indices(coords):
    """Axis = thinnest extent (empty/wedge direction); return the other two."""
    span = coords.max(axis=0) - coords.min(axis=0)
    axis = int(np.argmin(span))
    return [i for i in range(3) if i != axis]


def circle_center(ip):
    """Algebraic (Kasa) circle fit: x^2+y^2 = a*x + b*y + c -> centre (a/2, b/2)."""
    x, y = ip[:, 0], ip[:, 1]
    A = np.column_stack((x, y, np.ones_like(x)))
    rhs = x**2 + y**2
    a, b, _ = np.linalg.lstsq(A, rhs, rcond=None)[0]
    return np.array([0.5 * a, 0.5 * b])


def theta_from_upstream(coords, ip_idx, centre, flow_ip):
    """Angle (deg, 0..180) between each face's radial direction and the
    UPSTREAM direction (-flow). 0 = windward stagnation, 180 = leeward."""
    r = coords[:, ip_idx] - centre
    rn = np.linalg.norm(r, axis=1)
    rn[rn == 0.0] = 1.0
    rhat = r / rn[:, None]
    upstream = -flow_ip / np.linalg.norm(flow_ip)
    cosang = np.clip(rhat @ upstream, -1.0, 1.0)
    return np.degrees(np.arccos(cosang))


# --------------------------------- driver ---------------------------------- #
def main():
    q_inf = 0.5 * RHO_INF * U_INF ** 2
    print(f"[M20]  q_inf = {q_inf:.6g} Pa  (FOs already normalise Cp/Cf)")

    all_times = available_times()
    if not all_times:
        raise SystemExit(f"No time directories under {PP_DIR!r}")
    times = pick_times(all_times)
    if not times:
        raise SystemExit("No snapshots selected (check TIMES / TIME_RTOL).")

    # --- geometry established once (coordinates are identical across times) --
    flow3 = np.asarray(FLOW_DIR, dtype=float)
    geom = {}   # lazily filled: ip_idx, centre, flow_ip

    def ensure_geometry(coords):
        if geom:
            return
        ip_idx = in_plane_indices(coords)
        flow_ip = flow3[ip_idx]
        if np.linalg.norm(flow_ip) < 1e-12:
            raise SystemExit(
                "FLOW_DIR has no in-plane component (freestream parallel to the "
                "cylinder axis?). Set FLOW_DIR to the true freestream direction."
            )
        centre = circle_center(coords[:, ip_idx])
        geom.update(ip_idx=ip_idx, centre=centre, flow_ip=flow_ip)
        print(f"  in-plane = {ip_idx}, centre = "
              f"({centre[0]:.4g}, {centre[1]:.4g}), theta=0 from freestream "
              f"({', '.join(f'{c:g}' for c in flow3)})")

    # --- one figure per quantity, snapshots overlaid ----------------------- #
    for q in QUANTITIES:
        fig, ax = plt.subplots(figsize=(9, 6))
        cmap = plt.get_cmap("turbo")
        n = max(len(times), 1)
        plotted = 0

        for i, (t, tdir) in enumerate(times):
            path = find_raw_file(tdir, q["field"])
            if path is None:
                print(f"  [skip {q['key']}] no {q['field']} file in {tdir}")
                continue

            coords, val = load_surface(path)
            ensure_geometry(coords)
            th = theta_from_upstream(coords, geom["ip_idx"],
                                     geom["centre"], geom["flow_ip"])

            y = val * q["scale"]
            if q["take_abs"]:
                y = np.abs(y)
            if q["key"] == "heatflux":
                y = Q_SIGN * y

            order = np.argsort(th)
            ax.plot(th[order], y[order], color=cmap(i / n), lw=1.6,
                    label=f"nexaFoam : {q['tag']} ~ {t:.2e}s")
            plotted += 1

        if plotted == 0:
            print(f"  [skip {q['key']}] nothing plotted")
            plt.close(fig)
            continue

        ax.set_xlabel(r"$\theta$ from stagnation point  [deg]", fontweight="bold")
        ax.set_ylabel(q["ylabel"], fontweight="bold")
        ax.set_xlim(0.0, THETA_MAX)
        if q["ylim"] is not None:
            ax.set_ylim(*q["ylim"])
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=9, framealpha=0.9)
        fig.tight_layout()

        out = f"m20_{q['key']}_theta_time.png"
        fig.savefig(out, dpi=200)
        print(f"  wrote {out}")
        plt.close(fig)


if __name__ == "__main__":
    main()
