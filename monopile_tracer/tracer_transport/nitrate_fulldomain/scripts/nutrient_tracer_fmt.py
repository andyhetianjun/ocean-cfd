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
TRACER_NZ = 100


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
    vq = np.interp(zq, binned["depth"].values, binned["val"].values)
    vq = np.maximum.accumulate(vq)          # nitrate cannot decrease with depth
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
    return (np.nan_to_num(g[0]).astype(np.float32),
            np.nan_to_num(g[1]).astype(np.float32),
            np.nan_to_num(g[2]).astype(np.float32))


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
        Xd = np.clip(X3 - u*np.float32(DT_S), xmin, xmax)
        Yd = np.clip(Y3 - v*np.float32(DT_S), ymin, ymax)
        Zd_raw = Z3 - w*np.float32(DT_S)
        Zd = np.where(Zd_raw > zmax, 2*zmax - Zd_raw, Zd_raw)
        Zd = np.where(Zd < zmin, 2*zmin - Zd, Zd).astype(np.float32)
        C = _trilinear_interpolate(x, y, z_tr, C, Xd, Yd, Zd)
        for lst, s in zip(slabs, grab(C)):
            lst.append(s.astype(np.float32))
        del u, v, w, Xd, Yd, Zd
        if (n+1) % 20 == 0:
            print(f"    step {n+1}/{len(TIMESTEPS)}", flush=True)
    return slabs


XZ_ASPECT     = 4.5
XY_ASPECT     = 1.5   # y stretched, so the pile draws as an ellipse
CBAR_LABELPAD = 3
NUTRICLINE_BASE = 30.0   # set at runtime from the profile
NUTRICLINE_PEAK = 17.0   # set at runtime, depth of max |dC/dz|
CBAR_RANGES = {
    ("xz", None): (2.0,  8.0,  1.0),
    ("xy",  0):   (0.08, 0.24, 0.02),
    ("xy", 10):   (0.0,  1.2,  0.2),
    ("xy", 12):   (0.0,  2.0,  0.25),
    ("xy", 17):   (1.5,  5.0,  0.5),
    ("xy", 20):   (2.0,  6.0,  0.5),
    ("xy", 30):   (7.0, 10.0,  0.5),
    ("xy", 40):   (8.0, 10.0,  0.25),
}


def animate(mode, slice_val, slabs, x, y, z, filename, title_str="", fps=24, units_label=""):
    if mode == "xz":
        depth = Z_SURFACE - z                      # z=+20 surface -> depth 0; z=-20 -> depth 40
        H, V = np.meshgrid(x, depth)
        figsize, xl, yl = (14, 3.5), "x (m)", "z (m)"
        rkey = ("xz", None)
    else:
        H, V = np.meshgrid(x, y)
        figsize, xl, yl = ((14, 3.5) if FULL_DOMAIN else (9, 4)), "x (m)", "y (m)"
        _dep = Z_SURFACE - slice_val
        _cand = [k[1] for k in CBAR_RANGES if k[0] == "xy"]
        rkey = ("xy", min(_cand, key=lambda v: abs(v - _dep)))
    rng = CBAR_RANGES.get(rkey)
    if rng is not None:
        lo, hi, step = rng
    else:
        lo = float(min(s.min() for s in slabs)); hi = float(max(s.max() for s in slabs))
        if hi - lo < 1e-9:
            lo, hi = lo - 0.05*abs(lo) - 1e-6, hi + 0.05*abs(hi) + 1e-6
        raw = (hi - lo) / 6.0
        mag = 10.0 ** np.floor(np.log10(raw))
        step = next(m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag)
        lo = np.floor(lo / step) * step
        hi = np.ceil(hi / step) * step
    lo = max(lo, 0.0)
    cticks = np.round(np.arange(lo, hi + step * 0.5, step), 10)
    norm = Normalize(vmin=lo, vmax=hi)
    _tot = sum(s.size for s in slabs)
    _clip = sum(int((s < lo).sum()) + int((s > hi).sum()) for s in slabs)
    _ext = "both" if _clip else "neither"
    print(f"      colour range {lo:.3f} .. {hi:.3f}   clipped {100.0*_clip/_tot:.2f}%", flush=True)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.pcolormesh(H, V, slabs[0], cmap="viridis", norm=norm, shading="auto")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_title(title_str, fontsize=12, fontweight="bold", pad=8)
    cax = ax.inset_axes([1.02, 0, 0.015, 1], transform=ax.transAxes)
    _cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"),
                       cax=cax, ticks=cticks, extend=_ext)
    _cb.set_label(units_label, labelpad=CBAR_LABELPAD)
    fig.get_layout_engine().set(rect=(0, 0, 0.88, 1))
    ax.set_xlim(FULL_X[0], FULL_X[-1]) if FULL_DOMAIN else ax.set_xlim(x[0], x[-1])
    if mode == "xz":
        ax.set_ylim(40.0, 0.0)          # 0 (surface) top, 40 (bed) bottom
        ax.set_yticks(np.arange(0, 41, 5))
        ax.set_aspect(XZ_ASPECT, adjustable="box")
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
        ax.add_patch(patches.Circle((CYL_X, CYL_Y), CYL_R, facecolor="black",
                                    edgecolor="black", zorder=5))
        ax.set_aspect(XY_ASPECT, adjustable="box")  # y stretched for slide shape
    ax.autoscale(False)
    _tx, _ta = (0.98, "right") if mode == "xz" else (0.02, "left")
    tt = ax.text(_tx, 0.05, "", transform=ax.transAxes, ha=_ta, va="bottom",
                 color="white", fontsize=13,
                 bbox=dict(facecolor="black", alpha=0.5, edgecolor="none"))
    def update(i):
        im.set_array(slabs[i].ravel()); tt.set_text(f"t = {int(TIMESTEPS[0]) + int(i*DT_S)} s")
        return (im, tt)
    FuncAnimation(fig, update, frames=len(slabs), blit=False).save(
        filename, writer=PillowWriter(fps=fps))
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
    pdp, pvl = azmp_profile(AZMP_CSV, NUTRIENT)
    global NUTRICLINE_BASE, NUTRICLINE_PEAK
    _cp, _ = build_initial_field(z_tr, pdp, pvl, 1, 1)
    _c = _cp[:, 0, 0][::-1]                       # surface -> bed
    _d = (Z_SURFACE - z_tr)[::-1]
    _g = np.gradient(_c, _d)
    _pk = int(np.argmax(_g))                      # steepest point
    _thr = 0.2 * _g[_pk]
    _below = np.where(_g[_pk:] < _thr)[0]
    NUTRICLINE_BASE = float(_d[_pk + _below[0]]) if _below.size else float(_d[-1])
    NUTRICLINE_PEAK = float(_d[_pk])
    print(f"nutricline: peak gradient at {_d[_pk]:.1f} m, "
          f"base at {NUTRICLINE_BASE:.1f} m", flush=True)
    C0, depth_equiv = build_initial_field(z_tr, pdp, pvl, ny, nx)
    print(f"IC {NUTRIENT}: {C0.min():.3f} .. {C0.max():.3f}", flush=True)
    lab = "Nitrate (mmol m$^{-3}$)"; nut = "Nitrate"; fp = "nitrate"
    # NUTRIENT stays "nitrate" for the AZMP column lookup; only display/filenames change
    jobs = [("xz", 56.0, f"{fp}_xz_y56.gif",
             f"{nut}  |  vertical X-Z slice at y = 56 m (centre - 4 m, pile edge)"),
            ("xz", 58.0, f"{fp}_xz_y58.gif",
             f"{nut}  |  vertical X-Z slice at y = 58 m (centre - 2 m)"),
            ("xz", 60.0, f"{fp}_xz_y60.gif",
             f"{nut}  |  vertical X-Z slice at y = 60 m (centre, through monopile)"),
            ("xz", 62.0, f"{fp}_xz_y62.gif",
             f"{nut}  |  vertical X-Z slice at y = 62 m (centre + 2 m)"),
            ("xz", 64.0, f"{fp}_xz_y64.gif",
             f"{nut}  |  vertical X-Z slice at y = 64 m (centre + 4 m, pile edge)"),
            ("xy", float(z_tr[-1]), f"{fp}_xy_d00_surface.gif",
             f"{nut}  |  horizontal X-Y plane at z = 0 m (surface)"),
            ("xy", Z_SURFACE - 12.0, f"{fp}_xy_d12_nutriclinetop.gif",
             f"{nut}  |  horizontal X-Y plane at z = 12 m (top of nutricline)"),
            ("xy", Z_SURFACE - NUTRICLINE_PEAK,
             f"{fp}_xy_d{int(round(NUTRICLINE_PEAK)):02d}_nutriclinepeak.gif",
             f"{nut}  |  horizontal X-Y plane at z = {NUTRICLINE_PEAK:.0f} m "
             f"(steepest gradient)"),
            ("xy", Z_SURFACE - NUTRICLINE_BASE,
             f"{fp}_xy_d{int(round(NUTRICLINE_BASE)):02d}_nutriclinebase.gif",
             f"{nut}  |  horizontal X-Y plane at z = {NUTRICLINE_BASE:.0f} m "
             f"(base of nutricline)"),
            ("xy", 10.0, f"{fp}_xy_d10_subsurface.gif",
             f"{nut}  |  horizontal X-Y plane at z = 10 m (subsurface)"),
            ("xy", 0.0, f"{fp}_xy_d20_middle.gif",
             f"{nut}  |  horizontal X-Y plane at z = 20 m (middle)"),
            ("xy", -10.0, f"{fp}_xy_d30_subbottom.gif",
             f"{nut}  |  horizontal X-Y plane at z = 30 m (sub-bottom)"),
            ("xy", float(z_tr[0]), f"{fp}_xy_d40_bottom.gif",
             f"{nut}  |  horizontal X-Y plane at z = 40 m (bottom)")]
    import time
    CACHE = "data/slabs_cache"
    sig = "|".join([os.path.basename(VEL_DIR.rstrip("/")), str(TRACER_NZ), TIMESTEPS[0], TIMESTEPS[-1],
                    str(len(TIMESTEPS))] + [f"{m}:{sv}" for m, sv, _, _ in jobs])
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
