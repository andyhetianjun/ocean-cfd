#!/usr/bin/env python
"""
Nutrient tracer advection through the LES cylinder wake.

Physics core is Rhys' (semi-Lagrangian back-trace + trilinear interpolation +
AZMP profile IC). Changed for this case:
  - velocity from pre-gridded .npz per timestep (not ParaView CSVs)
  - z grid is the LES mesh's own 84 levels, so velocity needs no vertical interp
  - AZMP 3-95 m profile stretched onto the 40 m domain (Rhys' squash, rescaled)
  - adds an X-Y (plan) animator, which Rhys' code did not have
  - solves once, shared by all six views
  - colour scale per slice, from its own range across all frames
"""

import os
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation, PillowWriter

VEL_DIR    = "velocity_timeseries_z84"
AZMP_CSV   = "data/AZMP_Discrete_Occupations_Sections.csv"
NUTRIENT   = "nitrate"
UNITS      = {"nitrate": "Nitrate (mmol m$^{-3}$)",
              "phosphate": "Phosphate (mmol m$^{-3}$)",
              "silicate": "Silicate (mmol m$^{-3}$)"}

TIMESTEPS  = [str(t) for t in range(399, 459)]
DT_S       = 1.0

Z_SURFACE  = 20.0
Z_BOTTOM   = -20.0
AZMP_TOP   = 3.0
AZMP_BOT   = 95.0

CYL_X, CYL_Y, CYL_R = 46.0, 60.0, 4.0


def azmp_profile(csv_path, nutrient="nitrate"):
    """HL5, summer, binned at standard depths, spline-smoothed. Rhys' approach."""
    d = pd.read_csv(csv_path)
    d["date"] = pd.to_datetime(d["date"], format="%Y-%m-%d", errors="coerce")
    d["depth"] = pd.to_numeric(d["depth"], errors="coerce")
    m = (d["date"].notna() & d["date"].dt.month.isin([6, 7, 8, 9]) &
         (d["station"] == "HL5") & (d["depth"] <= 200))
    s = d.loc[m].copy()
    s[nutrient] = pd.to_numeric(s[nutrient], errors="coerce")

    rows = []
    for zz in [3, 10, 20, 30, 40, 50, 60, 80, 95]:
        b = s[(s["depth"] >= zz - 3) & (s["depth"] <= zz + 3)]
        rows.append({"depth": zz, "val": b[nutrient].mean(skipna=True)})
    binned = pd.DataFrame(rows).dropna()

    zq = np.linspace(AZMP_TOP, AZMP_BOT, 400)
    vq = UnivariateSpline(binned["depth"].values, binned["val"].values, s=0.0)(zq)
    return zq, np.clip(vq, 0, None)


def build_initial_field(z, prof_depth, prof_val, ny, nx):
    """
    The domain is 40 m but the summer nutricline sits at 40-100 m, so a direct
    depth mapping leaves almost no vertical gradient in the box. We stretch the
    AZMP 3-95 m profile across the 40 m domain, as Rhys did for his 6 m case.
    The vertical scale is therefore equivalent depth, not literal depth.
    """
    frac = (Z_SURFACE - z) / (Z_SURFACE - Z_BOTTOM)
    depth_equiv = AZMP_TOP + frac * (AZMP_BOT - AZMP_TOP)
    prof = np.interp(depth_equiv, prof_depth, prof_val)
    return (prof[:, None, None] * np.ones((1, ny, nx))).astype(np.float64), depth_equiv


def load_velocity(t):
    d = np.load(os.path.join(VEL_DIR, f"vel_{t}.npz"))
    g = d["grid_vel"]
    return (np.nan_to_num(g[0]).astype(np.float64),
            np.nan_to_num(g[1]).astype(np.float64),
            np.nan_to_num(g[2]).astype(np.float64))


def load_grid():
    d = np.load(os.path.join(VEL_DIR, f"vel_{TIMESTEPS[0]}.npz"))
    return d["x"].astype(float), d["y"].astype(float), d["z"].astype(float), d["counts"]


def _trilinear_interpolate(x, y, z, F, xq, yq, zq):
    xq = np.clip(xq, x[0], x[-1]); yq = np.clip(yq, y[0], y[-1]); zq = np.clip(zq, z[0], z[-1])
    nx, ny, nz = x.size, y.size, z.size
    ix = np.clip(np.searchsorted(x, xq, side="right") - 1, 0, nx - 2)
    iy = np.clip(np.searchsorted(y, yq, side="right") - 1, 0, ny - 2)
    iz = np.clip(np.searchsorted(z, zq, side="right") - 1, 0, nz - 2)
    x0, x1 = x[ix], x[ix + 1]; y0, y1 = y[iy], y[iy + 1]; z0, z1 = z[iz], z[iz + 1]
    tx = (xq - x0) / np.where(x1 > x0, x1 - x0, 1.0)
    ty = (yq - y0) / np.where(y1 > y0, y1 - y0, 1.0)
    tz = (zq - z0) / np.where(z1 > z0, z1 - z0, 1.0)
    shape = xq.shape
    ixf, iyf, izf = ix.ravel(), iy.ravel(), iz.ravel()
    G = lambda k, j, i: F[izf + k, iyf + j, ixf + i]
    txf, tyf, tzf = tx.ravel(), ty.ravel(), tz.ravel()
    c00 = G(0,0,0)*(1-txf) + G(0,0,1)*txf
    c01 = G(0,1,0)*(1-txf) + G(0,1,1)*txf
    c10 = G(1,0,0)*(1-txf) + G(1,0,1)*txf
    c11 = G(1,1,0)*(1-txf) + G(1,1,1)*txf
    c0 = c00*(1-tyf) + c01*tyf
    c1 = c10*(1-tyf) + c11*tyf
    return (c0*(1-tzf) + c1*tzf).reshape(shape)


def evolve(C0, x, y, z, jobs):
    """
    Run the advection once, keeping only the slices each view needs.

    Storing every full 3D frame costs len(TIMESTEPS) x nz*ny*nx*8 bytes, which
    is ~8 GB at a 800x240x84 grid. Extracting slices as we go keeps it in the
    tens of MB regardless of grid size.
    """
    C = C0.copy()
    nz, ny, nx = C.shape
    Z3, Y3, X3 = np.meshgrid(z, y, x, indexing="ij")
    xmin, xmax = x[0], x[-1]; ymin, ymax = y[0], y[-1]; zmin, zmax = z[0], z[-1]

    idx = []
    for mode, sv, fn in jobs:
        if mode == "xz":
            idx.append(("xz", int(np.argmin(np.abs(y - sv)))))
        else:
            idx.append(("xy", int(np.argmin(np.abs(z - sv)))))

    def grab(F):
        return [F[:, i, :].copy() if m == "xz" else F[i, :, :].copy() for m, i in idx]

    slabs = [[s] for s in grab(C)]
    for t in TIMESTEPS:
        u, v, w = load_velocity(t)
        Xd = np.clip(X3 - u*DT_S, xmin, xmax)
        Yd = np.clip(Y3 - v*DT_S, ymin, ymax)
        Zd = np.clip(Z3 - w*DT_S, zmin, zmax)
        C = _trilinear_interpolate(x, y, z, C, Xd, Yd, Zd)
        for lst, s in zip(slabs, grab(C)):
            lst.append(s)
    return slabs


def animate(mode, slice_val, slabs, x, y, z, filename, title_str="", fps=5, units_label=""):
    """
    Colour scale is taken from the min/max of THIS slice across ALL frames, so
    each panel uses its full colour range. A global scale washes the X-Y views
    out: mid-depth only spans ~4-5 of a 0-10 field, so everything lands in one
    shade. Consequence: scales differ between panels, so they cannot be
    compared to each other by eye.
    """
    if mode == "xz":
        H, V = np.meshgrid(x, z)
        figsize, xl, yl = (8, 4), "x (m)", "z (m)   [surface at top]"
    else:
        H, V = np.meshgrid(x, y)
        figsize, xl, yl = (9, 4), "x (m)", "y (m)"

    lo = float(min(s.min() for s in slabs))
    hi = float(max(s.max() for s in slabs))
    if hi - lo < 1e-9:
        lo, hi = lo - 0.05*abs(lo) - 1e-6, hi + 0.05*abs(hi) + 1e-6
    norm = Normalize(vmin=lo, vmax=hi)
    print(f"      colour range {lo:.3f} .. {hi:.3f}", flush=True)

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.pcolormesh(H, V, slabs[0], cmap="viridis", norm=norm, shading="auto")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_title(title_str, fontsize=12, fontweight="bold", pad=8)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"),
                 ax=ax, label=units_label, pad=0.005)
    # view limits must be set explicitly: pcolormesh sets DATA limits, but view
    # limits only autoscale at draw time, so freezing autoscale first leaves the
    # view at the default (0,1) and the plot renders blank.
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(V.min(), V.max())

    if mode == "xz" and abs(slice_val - CYL_Y) < CYL_R:
        half = np.sqrt(CYL_R**2 - (slice_val - CYL_Y)**2)
        ax.add_patch(patches.Rectangle((CYL_X - half, z[0]), 2*half, z[-1] - z[0],
                                       facecolor="black", edgecolor="black", zorder=5))
    elif mode == "xy":
        ax.add_patch(patches.Circle((CYL_X, CYL_Y), CYL_R, facecolor="black",
                                    edgecolor="white", zorder=5))
    ax.autoscale(False)
    tt = ax.text(0.98, 0.05, "", transform=ax.transAxes, ha="right", va="bottom",
                 color="white", fontsize=10,
                 bbox=dict(facecolor="black", alpha=0.5, edgecolor="none"))

    def update(i):
        im.set_array(slabs[i].ravel())
        tt.set_text(f"t = {int(i*DT_S)} s")
        return (im, tt)

    FuncAnimation(fig, update, frames=len(slabs), blit=False).save(
        filename, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return filename


def main():
    x, y, z, counts = load_grid()
    nz, ny, nx = len(z), len(y), len(x)
    print(f"grid {nx} x {ny} x {nz}   z {z[0]:.3f} .. {z[-1]:.3f}", flush=True)

    pdp, pvl = azmp_profile(AZMP_CSV, NUTRIENT)
    C0, depth_equiv = build_initial_field(z, pdp, pvl, ny, nx)
    print(f"IC {NUTRIENT}: {C0.min():.3f} .. {C0.max():.3f}  "
          f"(equiv depth {depth_equiv[-1]:.1f} m at z={z[-1]:.1f} "
          f"to {depth_equiv[0]:.1f} m at z={z[0]:.1f})", flush=True)

    lab = UNITS[NUTRIENT]
    nut = NUTRIENT.capitalize()
    jobs = [("xz", 60.0, f"{NUTRIENT}_xz_y60_centre.gif",
             f"{nut}  |  vertical X-Z slice at y = 60 m (centre, through monopile)"),
            ("xz", 52.0, f"{NUTRIENT}_xz_y52_minusD.gif",
             f"{nut}  |  vertical X-Z slice at y = 52 m (centre - D)"),
            ("xz", 68.0, f"{NUTRIENT}_xz_y68_plusD.gif",
             f"{nut}  |  vertical X-Z slice at y = 68 m (centre + D)"),
            ("xy", float(z[-1]), f"{NUTRIENT}_xy_top.gif",
             f"{nut}  |  horizontal X-Y plane at top Z (z = +20 m, surface)"),
            ("xy", float(z[nz//2]), f"{NUTRIENT}_xy_middle.gif",
             f"{nut}  |  horizontal X-Y plane at mid Z (z = 0 m)"),
            ("xy", float(z[0]), f"{NUTRIENT}_xy_bottom.gif",
             f"{nut}  |  horizontal X-Y plane at bottom Z (z = -20 m)")]

    import time
    t0 = time.time()
    print(f"advecting {len(TIMESTEPS)} steps (once, shared by all views)...", flush=True)
    all_slabs = evolve(C0, x, y, z, [(m, sv, fn) for m, sv, fn, _ in jobs])
    print(f"  done in {time.time()-t0:.0f}s", flush=True)

    for (mode, sv, fn, ttl), slabs in zip(jobs, all_slabs):
        t0 = time.time()
        print(f"  {fn}", flush=True)
        animate(mode, sv, slabs, x, y, z, fn, title_str=ttl, units_label=lab)
        print(f"      saved ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
