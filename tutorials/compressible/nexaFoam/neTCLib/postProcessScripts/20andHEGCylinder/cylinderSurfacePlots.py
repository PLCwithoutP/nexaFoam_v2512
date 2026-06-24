#!/usr/bin/env python3
"""
cylinder_surface_plots.py

Surface distributions vs. circumferential angle theta for the nexaFoam cylinder
validation cases (Mach-20 N2 half-cylinder and HEG air quarter-cylinder).

For each requested quantity it reads the corresponding raw 'surfaces' sample
written by the nexa* function objects on the `cylinder` patch and plots it
against theta measured from the stagnation point.

    theta(face) = angle of the face-centre about the cylinder axis,
                  referenced to the stagnation (max-pressure) face,
                  taken UNSIGNED so it runs 0 -> 180 deg
                  (a 90-deg quarter body simply tops out near 90).

The cylinder centre is recovered by an algebraic (Kasa) circle fit to the
in-plane face coordinates, so nothing about the centre/radius/orientation is
hard-coded. The cylinder axis is auto-detected as the thin (empty/wedge)
direction.

The FOs already normalise Cp, Cf and St, so those are plotted as-read. Only the
two dimensional fields (pressure, wall heat flux) get a unit scale per case.

Usage
    python cylinder_surface_plots.py HEG
    python cylinder_surface_plots.py M20
    (or set CASE below and run with no argument)
"""

import glob
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# =============================== SELECT CASE ================================ #
CASE = "HEG"            # "HEG" or "M20"; overridden by argv[1] if given
# =========================================================================== #


# ----------------------------- CASE CONFIG --------------------------------- #
# Each quantity:
#   key        : short id (used for the output png name)
#   field      : raw filename stem -> postProcessing/<FUNC>/<time>/<field>_<surf>.raw
#   scale      : multiply the sampled value (unit conversion)
#   take_abs   : abs() before plotting (heat flux sign convention)
#   ylabel/ylim/title
#   refs       : list of (csv_file, linestyle, marker, color, label)
#                CSV columns are  theta_deg, value   (plotted only if file exists)

COMMON_SURF = "cylinder"          # surface name inside the `surfaces` FO
FUNC_DIR    = "cylinderSurface"   # name of the `surfaces` function object

CASES = {
    "HEG": dict(
        title=r"reacting air  (HEG cylinder)",
        rho_inf=1.547e-3, U_inf=5956.0,
        theta_max=90.0,                       # quarter cylinder
        raw_glob_dir=f"postProcessing/{FUNC_DIR}/*",
        surf=COMMON_SURF,
        quantities=[
            dict(key="pressure", field="p",            scale=1e-3, take_abs=False,
                 ylabel=r"Surface pressure  [kPa]", ylim=(0.0, 55.0),
                 refs=[("su2nemo_p.csv", "-",   "None", "C1", "SU2-NEMO"),
                       ("exp_p.csv",     "none", "o",   "k",  "Experimental")]),
            dict(key="heatflux", field="wallHeatFlux", scale=1e-6, take_abs=True,
                 ylabel=r"Surface heat flux  [MW/m$^2$]", ylim=(0.0, 8.0),
                 refs=[("su2nemo_q.csv",  "-",   "None", "C1", "SU2-NEMO"),
                       ("exp_q.csv",      "none", "o",   "k",  "Experimental"),
                       ("nompelis_q.csv", "--",  "None", "C2", "Nompelis no Cat."),
                       ("lani_q.csv",     "--",  "None", "C0", "Lani no Cat.")]),
            dict(key="cp", field="cp_surface",    scale=1.0, take_abs=False,
                 ylabel=r"Pressure coefficient  $C_p$", ylim=(0.0, 2.0), refs=[]),
            dict(key="cf", field="cf_surface",    scale=1.0, take_abs=False,
                 ylabel=r"Skin-friction coefficient  $C_f$", ylim=(0.0, None), refs=[]),
            dict(key="st", field="StantonNumber", scale=1.0, take_abs=False,
                 ylabel=r"Stanton number  $St$", ylim=(0.0, None), refs=[]),
        ],
    ),
    "M20": dict(
        title=r"non-reacting $\mathrm{N_2}$  (Mach-20 cylinder)",
        rho_inf=1.363e-5, U_inf=6047.0,
        theta_max=180.0,                      # half cylinder
        raw_glob_dir=f"postProcessing/{FUNC_DIR}/*",
        surf=COMMON_SURF,
        quantities=[
            dict(key="pressure", field="p",            scale=1.0,  take_abs=False,
                 ylabel=r"Surface pressure  [Pa]", ylim=(0.0, None), refs=[]),
            dict(key="heatflux", field="wallHeatFlux", scale=1e-4, take_abs=True,
                 ylabel=r"Surface heat flux  [W/cm$^2$]", ylim=(0.0, 15.0),
                 refs=[("dsmcFoam_q.csv", "none", "+", "k", "DSMC: dsmcFoam"),
                       ("hy2Foam_q.csv",  "-",    "None", "k", "CFD: hy2Foam")]),
            dict(key="cp", field="cp_surface",    scale=1.0, take_abs=False,
                 ylabel=r"Pressure coefficient  $C_p$", ylim=(0.0, 2.0),
                 refs=[("dsmcFoam_cp.csv", "none", "+", "k", "DSMC: dsmcFoam"),
                       ("hy2Foam_cp.csv",  "-",    "None", "k", "CFD: hy2Foam")]),
            dict(key="cf", field="cf_surface",    scale=1.0, take_abs=False,
                 ylabel=r"Skin-friction coefficient  $C_f$", ylim=(0.0, 0.06),
                 refs=[("dsmcFoam_cf.csv", "none", "+", "k", "DSMC: dsmcFoam"),
                       ("hy2Foam_cf.csv",  "-",    "None", "k", "CFD: hy2Foam")]),
            dict(key="st", field="StantonNumber", scale=1.0, take_abs=False,
                 ylabel=r"Stanton number  $St$", ylim=(0.0, None), refs=[]),
        ],
    ),
}
# --------------------------------------------------------------------------- #


# ------------------------------ I/O helpers -------------------------------- #
def time_from_path(path):
    name = os.path.basename(os.path.dirname(path))
    try:
        return float(name)
    except ValueError:
        return float("nan")


def latest_raw(raw_glob_dir, field, surf):
    """Newest <field>_<surf>.raw under postProcessing/<func>/<time>/."""
    pattern = os.path.join(raw_glob_dir, f"{field}_{surf}.raw")
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


def theta_deg(coords, ip_idx, centre, theta_stag_deg):
    """Unsigned angular distance (deg, 0..180) from the stagnation direction."""
    ip = coords[:, ip_idx]
    ang = np.degrees(np.arctan2(ip[:, 1] - centre[1], ip[:, 0] - centre[0]))
    d = ang - theta_stag_deg
    d = (d + 180.0) % 360.0 - 180.0          # wrap to (-180, 180]
    return np.abs(d)


# --------------------------------- driver ---------------------------------- #
def main():
    case_name = (sys.argv[1] if len(sys.argv) > 1 else CASE).upper()
    if case_name not in CASES:
        raise SystemExit(f"Unknown case '{case_name}'. Choose from {list(CASES)}.")
    cfg = CASES[case_name]

    q_inf = 0.5 * cfg["rho_inf"] * cfg["U_inf"] ** 2
    print(f"[{case_name}]  q_inf = {q_inf:.6g} Pa  (FOs already normalise Cp/Cf/St)")

    # --- establish theta reference once, from the pressure sample ---------- #
    p_path = latest_raw(cfg["raw_glob_dir"], "p", cfg["surf"])
    pcoords, pvals = load_surface(p_path)
    ip_idx = in_plane_indices(pcoords)
    centre = circle_center(pcoords[:, ip_idx])
    ip = pcoords[:, ip_idx]
    stag_ang = np.degrees(np.arctan2(ip[:, 1][np.argmax(pvals)] - centre[1],
                                     ip[:, 0][np.argmax(pvals)] - centre[0]))
    print(f"  axis index = {[i for i in range(3) if i not in ip_idx][0]}, "
          f"in-plane = {ip_idx}, centre = ({centre[0]:.4g}, {centre[1]:.4g})")
    print(f"  stagnation face: p_max = {pvals.max():.6g} Pa, "
          f"theta_stag = {stag_ang:.2f} deg")
    th_chk = theta_deg(pcoords, ip_idx, centre, stag_ang)
    print(f"  theta span = [{th_chk.min():.1f}, {th_chk.max():.1f}] deg  "
          f"(expected ~0..{cfg['theta_max']:.0f}; if it stops near half, the patch "
          f"is symmetric about the stagnation line and theta is folded)")

    # --- one figure per requested quantity --------------------------------- #
    for q in cfg["quantities"]:
        try:
            raw = latest_raw(cfg["raw_glob_dir"], q["field"], cfg["surf"])
        except FileNotFoundError as e:
            print(f"  [skip {q['key']}] {e}")
            continue

        coords, val = load_surface(raw)
        th = theta_deg(coords, ip_idx, centre, stag_ang)

        y = val * q["scale"]
        if q["take_abs"]:
            y = np.abs(y)

        order = np.argsort(th)
        th, y = th[order], y[order]
        print(f"  {q['key']:9s} <- {os.path.basename(raw)}  "
              f"(N={val.size}, raw range [{val.min():.4g}, {val.max():.4g}])")

        fig, ax = plt.subplots(figsize=(7.0, 4.5))
        ax.plot(th, y, "-", color="r", lw=1.4, label="CFD: nexaFoam")

        for ref_file, ls, marker, color, label in q["refs"]:
            if ref_file and os.path.exists(ref_file):
                ref = np.loadtxt(ref_file, delimiter=",")
                if ref.ndim == 1:
                    ref = ref[None, :]
                ax.plot(ref[:, 0], ref[:, 1], linestyle=ls, marker=marker,
                        color=color, ms=5, lw=1.4, label=label)

        ax.set_xlabel(r"$\theta$ from stagnation point  [deg]")
        ax.set_ylabel(q["ylabel"])
        ax.set_xlim(0.0, cfg["theta_max"])
        lo, hi = q["ylim"]
        ax.set_ylim(bottom=lo if lo is not None else None,
                    top=hi if hi is not None else None)
        ax.set_title(cfg["title"])
        ax.legend(frameon=True)
        fig.tight_layout()

        out = f"{case_name.lower()}_{q['key']}_theta.png"
        fig.savefig(out, dpi=150)
        print(f"            -> wrote {out}")
        plt.close(fig)


if __name__ == "__main__":
    main()
