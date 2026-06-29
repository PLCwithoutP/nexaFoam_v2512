#!/usr/bin/env python3
"""
cylinder_surface_plots_m20.py

Surface distributions vs. circumferential angle theta for the nexaFoam
Mach-20 N2 half-cylinder validation case.

Three quantities are plotted, each overlaying three data sets:
    * pressure coefficient   Cp
    * skin-friction coeff.   Cf
    * surface heat flux      [W/cm^2]

Data sets / styles:
    nexaFoam   -> red solid line          (this work)
    dsmcFoam   -> black markers           (DSMC reference)
    hy2Foam    -> dark-blue dashed line   (CFD reference)

theta is measured from the windward stagnation point (0 deg) to the leeward
side (180 deg). The cylinder centre is recovered by an algebraic (Kasa) circle
fit to the face coordinates and theta is built purely from geometry + the known
freestream direction, so it is immune to spurious peaks in the solution (e.g. a
corner blow-up that would hijack a max-pressure anchor).

Usage
    python cylinder_surface_plots_m20.py
"""

import glob
import os
import numpy as np
import matplotlib.pyplot as plt

# =============================== CONFIG ==================================== #
RHO_INF  = 1.363e-5        # kg/m^3   (Table 2)
U_INF    = 6047.0          # m/s
FLOW_DIR = (1.0, 0.0, 0.0)  # freestream velocity direction in the mesh frame
THETA_MAX = 180.0          # half cylinder

FUNC_DIR = "coneSurface"   # name of the `surfaces` function object
                               # -> postProcessing/<FUNC_DIR>/<time>/<field>_<surf>.raw
SURF     = "cylinder"          # surface name inside the `surfaces` FO

RAW_GLOB_DIR = f"postProcessing/{FUNC_DIR}/*"

# Plot styles for the nexaFoam (this-work) curve.
NEXA_STYLE = dict(linestyle="-", marker="None", color="red", lw=1.4,
                  label="CFD: nexaFoam")

# Each quantity:
#   key      : short id (used for the output png name)
#   field    : raw filename stem -> <field>_<surf>.raw
#   scale    : multiply the sampled value (unit conversion)
#   take_abs : abs() before plotting (heat-flux sign convention)
#   refs     : list of (csv_file, linestyle, marker, color, label)
#              CSV columns are  theta_deg, value   (plotted only if file exists)
DSMC_STYLE = ("none", "o", "k")          # (linestyle, marker, color)
HY2_STYLE  = ("--",   "None", "darkblue")

QUANTITIES = [
    dict(key="cp", field="cp_surface", scale=1.0, take_abs=False,
         ylabel=r"Pressure coefficient  $C_p$",
         refs=[("dsmcFoam_cp.csv", *DSMC_STYLE, "DSMC: dsmcFoam"),
               ("hy2Foam_cp.csv",  *HY2_STYLE,  "CFD: hy2Foam")]),
    dict(key="cf", field="cf_surface", scale=1.0, take_abs=False,
         ylabel=r"Skin-friction coefficient  $C_f$",
         refs=[("dsmcFoam_cf.csv", *DSMC_STYLE, "DSMC: dsmcFoam"),
               ("hy2Foam_cf.csv",  *HY2_STYLE,  "CFD: hy2Foam")]),
    dict(key="heatflux", field="wallHeatFlux", scale=1e-4, take_abs=True,
         ylabel=r"Surface heat flux  [W/cm$^2$]",
         refs=[("dsmcFoam_q.csv", *DSMC_STYLE, "DSMC: dsmcFoam"),
               ("hy2Foam_q.csv",  *HY2_STYLE,  "CFD: hy2Foam")]),
]

TITLE = r"Reacting $\mathrm{N_2}$  (Mach-20 cylinder)"
# =========================================================================== #


# ------------------------------ I/O helpers -------------------------------- #
def time_from_path(path):
    name = os.path.basename(os.path.dirname(path))
    try:
        return float(name)
    except ValueError:
        return float("nan")


def latest_raw(field, surf):
    """Newest <field>_<surf>.raw under postProcessing/<func>/<time>/."""
    pattern = os.path.join(RAW_GLOB_DIR, f"{field}_{surf}.raw")
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
    r = coords[:, ip_idx] - centre                  # N x 2 (face - centre)
    rn = np.linalg.norm(r, axis=1)
    rn[rn == 0.0] = 1.0
    rhat = r / rn[:, None]
    upstream = -flow_ip / np.linalg.norm(flow_ip)   # points into the flow
    cosang = np.clip(rhat @ upstream, -1.0, 1.0)
    return np.degrees(np.arccos(cosang))


# --------------------------------- driver ---------------------------------- #
def main():
    q_inf = 0.5 * RHO_INF * U_INF ** 2
    print(f"[M20]  q_inf = {q_inf:.6g} Pa  (FOs already normalise Cp/Cf/St)")

    # --- establish geometry once, from the pressure sample's coordinates --- #
    # (coordinates never corrupt; only the *values* might, so theta is built
    #  from geometry + freestream direction, not from the pressure field.)
    pcoords, pvals = load_surface(latest_raw("p", SURF))
    ip_idx = in_plane_indices(pcoords)
    centre = circle_center(pcoords[:, ip_idx])

    flow3 = np.asarray(FLOW_DIR, dtype=float)
    flow_ip = flow3[ip_idx]
    if np.linalg.norm(flow_ip) < 1e-12:
        raise SystemExit(
            "FLOW_DIR has no in-plane component (freestream parallel to the "
            "cylinder axis?). Set FLOW_DIR to the true freestream direction."
        )

    axis_idx = [i for i in range(3) if i not in ip_idx][0]
    print(f"  axis index = {axis_idx}, in-plane = {ip_idx}, "
          f"centre = ({centre[0]:.4g}, {centre[1]:.4g})")
    print(f"  theta=0 fixed by freestream direction "
          f"({', '.join(f'{c:g}' for c in flow3)}) (geometric, not pressure-based)")

    # diagnostic: where does the pressure peak actually sit?
    th_p = theta_from_upstream(pcoords, ip_idx, centre, flow_ip)
    i_pmax = int(np.argmax(pvals))
    print(f"  theta span = [{th_p.min():.1f}, {th_p.max():.1f}] deg "
          f"(expected ~0..{THETA_MAX:.0f})")
    print(f"  max pressure {pvals.max():.4g} Pa sits at theta = {th_p[i_pmax]:.1f} deg"
          f"  -> {'OK (near stagnation)' if th_p[i_pmax] < 10 else 'WARNING: spurious peak away from stagnation (corner/blow-up?)'}")

    # --- one figure per requested quantity --------------------------------- #
    for q in QUANTITIES:
        try:
            raw = latest_raw(q["field"], SURF)
        except FileNotFoundError as e:
            print(f"  [skip {q['key']}] {e}")
            continue

        coords, val = load_surface(raw)
        th = theta_from_upstream(coords, ip_idx, centre, flow_ip)

        y = val * q["scale"]
        if q["take_abs"]:
            y = np.abs(y)

        order = np.argsort(th)
        th, y = th[order], y[order]
        print(f"  {q['key']:9s} <- {os.path.basename(raw)}  "
              f"(N={val.size}, raw range [{val.min():.4g}, {val.max():.4g}])")

        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        ax.plot(th, y, **NEXA_STYLE)

        for ref_file, ls, marker, color, label in q["refs"]:
            if ref_file and os.path.exists(ref_file):
                ref = np.loadtxt(ref_file, delimiter=",")
                if ref.ndim == 1:
                    ref = ref[None, :]
                ax.plot(ref[:, 0], ref[:, 1], linestyle=ls, marker=marker,
                        color=color, ms=5, lw=1.4, label=label)

        ax.set_xlabel(r"$\theta$ from stagnation point  [deg]")
        ax.set_ylabel(q["ylabel"])
        ax.set_xlim(0.0, THETA_MAX)
        # y-axis follows the data (CFD curve + any references).
        ax.relim()
        ax.autoscale(axis="y", tight=False)
        ax.margins(y=0.05)
        ax.set_title(TITLE)
        ax.legend(frameon=True)
        fig.tight_layout()

        out = f"m20_{q['key']}_theta.png"
        fig.savefig(out, dpi=150)
        print(f"            -> wrote {out}")
        plt.close(fig)


if __name__ == "__main__":
    main()
