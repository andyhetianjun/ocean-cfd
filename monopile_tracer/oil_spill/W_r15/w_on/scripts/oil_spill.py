#!/usr/bin/env python
"""
Fine in all three: 0.25 m x-y (velocity_timeseries_fine) + 400 z-levels (0.1 m,
interpolated from the 84 computed levels) + 197 consecutive timesteps (~2.5
shedding cycles). Vertical detail is interpolated, not computed. float32 in the
heavy path to keep memory bounded on a shared box.
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

VEL_DIR    = "data/velocity_full_025"
AZMP_CSV   = "data/AZMP_Discrete_Occupations_Sections.csv"
NUTRIENT   = "nitrate"
UNITS      = {"nitrate": "Nitrate (mmol m$^{-3}$)",
              "phosphate": "Phosphate (mmol m$^{-3}$)",
              "silicate": "Silicate (mmol m$^{-3}$)"}

TIMESTEPS  = [str(t) for t in range(1049, 1401)]
DT_S       = 1.0

Z_SURFACE  = 20.0
Z_BOTTOM   = -20.0
AZMP_TOP   = 3.0
AZMP_BOT   = 95.0

CYL_X, CYL_Y, CYL_R = 40.0, 60.0, 4.0
TRACER_NZ = 150

# --- oil spill config -------------------------------------------------
OIL_MASS_KG   = 10.0        # 0.0118 m3 at 850 kg/m3, 16.6 um over r=15 m
OIL_RADIUS    = 15.0        # m
OIL_OFFSET    = 24.0        # m, spill centre to pile centre
BEARING       = "W"         # N NE E SE S SW W NW
VERTICAL_MIX  = True        # case 2
_BEAR = {"N": (0, 1), "NE": (0.7071, 0.7071), "E": (1, 0), "SE": (0.7071, -0.7071),
         "S": (0, -1), "SW": (-0.7071, -0.7071), "W": (-1, 0), "NW": (-0.7071, 0.7071)}
OIL_CX = CYL_X + _BEAR[BEARING][0] * OIL_OFFSET
OIL_CY = CYL_Y + _BEAR[BEARING][1] * OIL_OFFSET
# ----------------------------------------------------------------------


def azmp_profile(csv_path, nutrient="nitrate"):
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


FULL_DOMAIN = True          # pad computed 20-220 window out to full 0-600 for plotting
FULL_X = (0.0, 600.0)       # full LES domain extent in x
FULL_Y = (0.0, 120.0)       # full LES domain extent in y

def build_initial_field(z, prof_depth, prof_val, ny, nx):
    frac = (Z_SURFACE - z) / (Z_SURFACE - Z_BOTTOM)
    depth_equiv = AZMP_TOP + frac * (AZMP_BOT - AZMP_TOP)
    prof = np.interp(depth_equiv, prof_depth, prof_val)
    return (prof[:, None, None] * np.ones((1, ny, nx))).astype(np.float32), depth_equiv


def load_velocity(t):
    d = np.load(os.path.join(VEL_DIR, f"vel_{t}.npz"))
    g = d["grid_vel"]
    w = np.nan_to_num(g[2]).astype(np.float32)
    if not VERTICAL_MIX:
        w = np.zeros_like(w)
    return (np.nan_to_num(g[0]).astype(np.float32),
            np.nan_to_num(g[1]).astype(np.float32),
            w)


def build_oil_field(x, y, z_tr):
    """Disc of oil in the surface cell. Returns (C in g/m3, diagnostics)."""
    nz, ny, nx = len(z_tr), len(y), len(x)
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0]); dz = float(z_tr[1] - z_tr[0])
    X, Y = np.meshgrid(x, y)
    inside = ((X - OIL_CX)**2 + (Y - OIL_CY)**2) <= OIL_RADIUS**2
    disc_area = np.pi * OIL_RADIUS**2
    conc = (OIL_MASS_KG * 1000.0) / (disc_area * dz)      # g/m3
    C = np.zeros((nz, ny, nx), dtype=np.float32)
    C[-1][inside] = conc                                   # z_tr[-1] is the surface
    grid_area = float(inside.sum()) * dx * dy
    return C, {"conc_gm3": conc, "areal_gm2": conc * dz,
               "frac_in_grid": grid_area / disc_area,
               "mass_in_grid_kg": conc * grid_area * dz / 1000.0}


def load_grid():
    d = np.load(os.path.join(VEL_DIR, f"vel_{TIMESTEPS[0]}.npz"))
    x = d["x"].astype(float); y = d["y"].astype(float)
    z_vel = d["z"].astype(float)
    if TRACER_NZ and TRACER_NZ != len(z_vel):
        z_tr = np.linspace(z_vel[0], z_vel[-1], TRACER_NZ)
    else:
        z_tr = z_vel
    return x, y, z_vel, z_tr, d["counts"]


def velocity_on_tracer_z(u, v, w, z_vel, z_tr):
    if z_tr.shape == z_vel.shape and np.allclose(z_tr, z_vel):
        return u, v, w
    iz = np.clip(np.searchsorted(z_vel, z_tr, side="right") - 1, 0, len(z_vel) - 2)
    z0 = z_vel[iz]; z1 = z_vel[iz + 1]
    t = ((z_tr - z0) / np.where(z1 > z0, z1 - z0, 1.0))[:, None, None].astype(np.float32)
    lift = lambda F: F[iz] * (1 - t) + F[iz + 1] * t
    return lift(u), lift(v), lift(w)


def _trilinear_interpolate(x, y, z, F, xq, yq, zq):
    xq = np.clip(xq, x[0], x[-1]); yq = np.clip(yq, y[0], y[-1]); zq = np.clip(zq, z[0], z[-1])
    nx, ny, nz = x.size, y.size, z.size
    ix = np.clip(np.searchsorted(x, xq, side="right") - 1, 0, nx - 2)
    iy = np.clip(np.searchsorted(y, yq, side="right") - 1, 0, ny - 2)
    iz = np.clip(np.searchsorted(z, zq, side="right") - 1, 0, nz - 2)
    x0, x1 = x[ix], x[ix + 1]; y0, y1 = y[iy], y[iy + 1]; z0, z1 = z[iz], z[iz + 1]
    tx = ((xq - x0) / np.where(x1 > x0, x1 - x0, 1.0)).astype(np.float32)
    ty = ((yq - y0) / np.where(y1 > y0, y1 - y0, 1.0)).astype(np.float32)
    tz = ((zq - z0) / np.where(z1 > z0, z1 - z0, 1.0)).astype(np.float32)
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


def evolve(C0, x, y, z_vel, z_tr, jobs):
    C = C0.copy()
    Z3, Y3, X3 = (a.astype(np.float32) for a in np.meshgrid(z_tr, y, x, indexing="ij"))
    xmin, xmax = x[0], x[-1]; ymin, ymax = y[0], y[-1]; zmin, zmax = z_tr[0], z_tr[-1]
    idx = []
    for mode, sv, fn in jobs:
        if mode == "xz":
            idx.append(("xz", int(np.argmin(np.abs(y - sv)))))
        else:
            idx.append(("xy", int(np.argmin(np.abs(z_tr - sv)))))
    def grab(F):
        return [F[:, i, :].copy() if m == "xz" else F[i, :, :].copy() for m, i in idx]
    slabs = [[s] for s in grab(C)]
    for n, t in enumerate(TIMESTEPS):
        u, v, w = load_velocity(t)
        u, v, w = velocity_on_tracer_z(u, v, w, z_vel, z_tr)
        Xd_raw = X3 - u*np.float32(DT_S)
        outside = Xd_raw < xmin
        Xd = np.clip(Xd_raw, xmin, xmax)
        Yd = np.clip(Y3 - v*np.float32(DT_S), ymin, ymax)
        Zd_raw = Z3 - w*np.float32(DT_S)
        Zd = np.where(Zd_raw > zmax, 2*zmax - Zd_raw, Zd_raw)
        Zd = np.where(Zd < zmin, 2*zmin - Zd, Zd).astype(np.float32)
        C = _trilinear_interpolate(x, y, z_tr, C, Xd, Yd, Zd)
        C[outside] = 0.0
        for lst, s in zip(slabs, grab(C)):
            lst.append(s.astype(np.float32))
        del u, v, w, Xd, Yd, Zd
        if (n+1) % 20 == 0:
            _dv = float(x[1]-x[0]) * float(y[1]-y[0]) * float(z_tr[1]-z_tr[0])
            print(f"    step {n+1}/{len(TIMESTEPS)}  mass={C.sum()*_dv/1000.0:.4f} kg",
                  flush=True)
    return slabs


def animate(mode, slice_val, slabs, x, y, z, filename, title_str="", fps=12, units_label=""):
    if mode == "xz":
        depth = Z_SURFACE - z                      # z=+20 surface -> depth 0; z=-20 -> depth 40
        H, V = np.meshgrid(x, depth)
        figsize, xl, yl = (8, 4), "x (m)", "z (m)"
    else:
        H, V = np.meshgrid(x, y)
        figsize, xl, yl = ((14, 3.5) if FULL_DOMAIN else (9, 4)), "x (m)", "y (m)"
    lo = float(min(s.min() for s in slabs)); hi = float(max(s.max() for s in slabs))
    lo = max(lo, 0.0)   # oil concentration cannot be negative
    if hi - lo < 1e-9:
        lo, hi = lo - 0.05*abs(lo) - 1e-6, hi + 0.05*abs(hi) + 1e-6
    raw = (hi - lo) / 6.0
    mag = 10.0 ** np.floor(np.log10(raw))
    step = next(m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag)
    lo = np.floor(lo / step) * step
    hi = np.ceil(hi / step) * step
    cticks = np.round(np.arange(lo, hi + step * 0.5, step), 10)
    norm = Normalize(vmin=lo, vmax=hi)
    print(f"      colour range {lo:.3f} .. {hi:.3f}", flush=True)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.pcolormesh(H, V, slabs[0], cmap="viridis", norm=norm, shading="auto")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_title(title_str, fontsize=12, fontweight="bold", pad=8)
    cax = ax.inset_axes([1.02, 0, 0.015, 1], transform=ax.transAxes)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"),
                 cax=cax, label=units_label, ticks=cticks)
    fig.get_layout_engine().set(rect=(0, 0, 0.88, 1))
    ax.set_xlim(FULL_X[0], FULL_X[-1]) if FULL_DOMAIN else ax.set_xlim(x[0], x[-1])
    if mode == "xz":
        ax.set_ylim(40.0, 0.0)          # 0 (surface) top, 40 (bed) bottom
        ax.set_yticks(np.arange(0, 41, 5))
    else:
        if FULL_DOMAIN:
            ax.set_ylim(FULL_Y[0], FULL_Y[-1])
        else:
            ax.set_ylim(V.min(), V.max())
    if mode == "xz" and abs(slice_val - CYL_Y) < CYL_R:
        half = np.sqrt(CYL_R**2 - (slice_val - CYL_Y)**2)
        ax.add_patch(patches.Rectangle((CYL_X - half, V.min()), 2*half, V.max() - V.min(),
                                       facecolor="black", edgecolor="black", zorder=5))
    elif mode == "xy":
        ax.add_patch(patches.Circle((CYL_X, CYL_Y), CYL_R * 1.1, facecolor="black",
                                    edgecolor="black", zorder=5))
        ax.set_aspect(1.5, adjustable="box")   # matched to nitrate plan views
    ax.autoscale(False)
    _tx, _ta = (0.98, "right") if mode == "xz" else (0.02, "left")
    tt = ax.text(_tx, 0.05, "", transform=ax.transAxes, ha=_ta, va="bottom",
                 color="white", fontsize=10,
                 bbox=dict(facecolor="black", alpha=0.5, edgecolor="none"))
    def update(i):
        im.set_array(slabs[i].ravel()); tt.set_text(f"t = {int(TIMESTEPS[0]) + int(i*DT_S)} s")
        return (im, tt)
    FuncAnimation(fig, update, frames=len(slabs), blit=False).save(
        os.path.join("results", filename), writer=PillowWriter(fps=fps))
    plt.close(fig)
    return filename


def pad_to_full(slabs, x, y, C0, mode, slice_idx):
    """Pad each computed frame out to the full LES domain, filling outside
    the computed window with the t=0 initial profile value. x is padded for
    all modes; y is padded only for xy plan views (where the vertical axis
    is y). Returns (padded_slabs, x_full, y_full)."""
    dx = x[1] - x[0]
    xf0 = FULL_X[0] + dx/2.0
    nx_full = int(round((FULL_X[1] - FULL_X[0]) / dx))
    x_full = xf0 + np.arange(nx_full) * dx
    j0 = int(round((x[0] - x_full[0]) / dx))          # computed window start in x
    if mode == "xz":
        # rows are depth; fill each row with its initial value; no y padding
        fill_col = C0[:, 0, 0][:, None]               # (nz,1)
        out = []
        for s in slabs:
            base = np.repeat(fill_col, nx_full, axis=1).astype(np.float32)
            base[:, j0:j0 + s.shape[1]] = s
            out.append(base)
        return out, x_full, y
    # xy plan view: single plane value, pad both x and y
    plane_val = float(C0[slice_idx, 0, 0])
    dy = y[1] - y[0]
    yf0 = FULL_Y[0] + dy/2.0
    ny_full = int(round((FULL_Y[1] - FULL_Y[0]) / dy))
    y_full = yf0 + np.arange(ny_full) * dy
    i0 = int(round((y[0] - y_full[0]) / dy))          # computed window start in y
    out = []
    for s in slabs:
        base = np.full((ny_full, nx_full), plane_val, dtype=np.float32)
        base[i0:i0 + s.shape[0], j0:j0 + s.shape[1]] = s
        out.append(base)
    return out, x_full, y_full


def main():
    x, y, z_vel, z_tr, counts = load_grid()
    ny, nx = len(y), len(x)
    print(f"grid {nx} x {ny}   velocity z: {len(z_vel)} levels   "
          f"tracer z: {len(z_tr)} levels ({(z_tr[-1]-z_tr[0])/(len(z_tr)-1):.3f} m)", flush=True)
    C0, oil = build_oil_field(x, y, z_tr)
    print(f"oil: r={OIL_RADIUS} m at ({OIL_CX:.1f}, {OIL_CY:.1f}) bearing {BEARING}", flush=True)
    print(f"     {oil['conc_gm3']:.1f} g/m3 in surface cell = {oil['areal_gm2']:.2f} g/m2", flush=True)
    print(f"     {oil['frac_in_grid']*100:.1f}% of disc inside grid, "
          f"{oil['mass_in_grid_kg']:.2f} of {OIL_MASS_KG} kg", flush=True)
    print(f"     vertical mixing: {'ON (case 2)' if VERTICAL_MIX else 'OFF (case 1, w=0)'}", flush=True)
    lab = "Oil (g m$^{-2}$)" if not VERTICAL_MIX else "Oil (g m$^{-3}$)"
    nut = "Oil"; fp = f"oil_{BEARING}_case2"
    # NUTRIENT stays "nitrate" for the AZMP column lookup; only display/filenames change
    jobs = [("xy", float(z_tr[-1]), f"{fp}_xy_z0_surface.gif",
             f"{nut}  |  X-Y plane at z = 0 m (surface), spill {BEARING} of monopile, with no buoyancy")]
    import time
    CACHE = "data/slabs_cache"
    import hashlib
    ich = hashlib.md5(C0.tobytes()).hexdigest()[:12]
    sig = "|".join([os.path.basename(VEL_DIR.rstrip("/")), str(TRACER_NZ), TIMESTEPS[0], TIMESTEPS[-1],
                    str(len(TIMESTEPS)), ich, str(VERTICAL_MIX)]
                   + [f"{m}:{sv}" for m, sv, _, _ in jobs])
    sigf = os.path.join(CACHE, "sig.txt")
    all_slabs = None
    if os.path.exists(sigf) and open(sigf).read() == sig:
        all_slabs = [list(np.load(os.path.join(CACHE, f"arr_{k}.npy")))
                     for k in range(len(jobs))]
        print(f"  loaded cached slabs from {CACHE}/ (advection skipped)", flush=True)
    elif os.path.exists(sigf):
        print(f"  {CACHE}/ exists but config changed, re-advecting", flush=True)
    if all_slabs is None:
        t0 = time.time()
        print(f"advecting {len(TIMESTEPS)} steps on fine tracer grid...", flush=True)
        all_slabs = evolve(C0, x, y, z_vel, z_tr, [(m, sv, fn) for m, sv, fn, _ in jobs])
        print(f"  done in {time.time()-t0:.0f}s", flush=True)
        t0 = time.time()
        os.makedirs(CACHE, exist_ok=True)
        for k, s in enumerate(all_slabs):
            np.save(os.path.join(CACHE, f"arr_{k}.npy"), np.stack(s))
        open(sigf, "w").write(sig)
        print(f"  cached slabs to {CACHE}/ ({time.time()-t0:.0f}s)", flush=True)
    if not VERTICAL_MIX:
        # case 1: all oil in one layer, so g/m3 * dz is exact areal density.
        # case 2 spreads vertically and needs a column integral, not a slice.
        _dz = float(z_tr[1] - z_tr[0])
        all_slabs = [[s * _dz for s in sl] for sl in all_slabs]
    for (mode, sv, fn, ttl), slabs in zip(jobs, all_slabs):
        t0 = time.time()
        print(f"  {fn}", flush=True)
        xplot = x; yplot = y
        if FULL_DOMAIN:
            siz = int(np.argmin(np.abs(z_tr - sv))) if mode == "xy" else 0
            slabs, xplot, yplot = pad_to_full(slabs, x, y, C0, mode, siz)
        animate(mode, sv, slabs, xplot, yplot, z_tr, fn, title_str=ttl, units_label=lab)
        print(f"      saved ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
