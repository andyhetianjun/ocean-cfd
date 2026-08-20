#!/usr/bin/env python

import numpy as np
import pandas as pd

import sys
from typing import Tuple, Dict, Optional
import math
import io
import os

from scipy.interpolate import UnivariateSpline

from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import animation
from matplotlib.colors import Normalize
from matplotlib.animation import PillowWriter


# ---------------------------------------------------
# Load Data
# ---------------------------------------------------

# Load nutrient data
nutrient_data = pd.read_csv("data/AZMP_Discrete_Occupations_Sections.csv")

# ---------------------------------------------------
# Nutrient Data Creation
# ---------------------------------------------------

# 1) Ensure datetime
nutrient_data['date'] = pd.to_datetime(nutrient_data['date'], format='%Y-%m-%d', errors='coerce')

# 2) Ensure numeric depth (so comparisons don't fail)
nutrient_data['depth'] = pd.to_numeric(nutrient_data['depth'], errors='coerce')

# 3) Summer Months
summer_months = [6, 7, 8, 9]  # Jun–Sep;

mask = (
    nutrient_data['date'].notna() &
    nutrient_data['date'].dt.month.isin(summer_months) &
    (nutrient_data['station'] == 'HL5') &
    (nutrient_data['depth'] <= 200)
)

nutrient_subset = nutrient_data.loc[mask].copy()

# 4) Bin means
def compute_binned_means(
    df,
    standard_depths=[3,10,20,30,40,50,60,80,95],
    bin_half_width=3,
    depth_col="depth",
    nitrate_col="nitrate",
    phosphate_col="phosphate",
    silicate_col="silicate"
):
    data = df.copy()
    for c in [depth_col, nitrate_col, phosphate_col, silicate_col]:
        data[c] = pd.to_numeric(data[c], errors="coerce")

    rows = []
    for z in standard_depths:
        mask = (data[depth_col] >= z - bin_half_width) & (data[depth_col] <= z + bin_half_width)
        sub = data.loc[mask]

        rows.append({
            "depth": z,
            "nitrate": sub[nitrate_col].mean(skipna=True),
            "phosphate": sub[phosphate_col].mean(skipna=True),
            "silicate": sub[silicate_col].mean(skipna=True),
            "n": int(len(sub))
        })

    return pd.DataFrame(rows)

nutrient_binned = compute_binned_means(nutrient_subset)

# 5) Interpolate
def interpolate_profiles(
    nutrient_binned: pd.DataFrame,
    z_min: float = 0,
    z_max: float = 100,
    n_points: int = 400,
    smoothing: float = 0.0,         # used by spline/splrep (e.g., 0–5)
    lowess_frac: float = 0.3        # used by lowess
) -> pd.DataFrame:
    z = np.linspace(z_min, z_max, n_points)
    out = {"depth": z}

    for var in ["nitrate", "phosphate", "silicate"]:
        s = nutrient_binned[var]
        valid = s.notna()
        if valid.sum() == 0:
            out[var] = np.full_like(z, np.nan, dtype=float)
            continue
        if valid.sum() == 1:
            out[var] = np.full_like(z, s.loc[valid].iloc[0], dtype=float)
            continue

        xp = nutrient_binned.loc[valid, "depth"].values
        fp = s.loc[valid].values
        order = np.argsort(xp)
        xp, fp = xp[order], fp[order]

        y = UnivariateSpline(xp, fp, s=smoothing)(z)

        out[var] = y

    return pd.DataFrame(out, columns=["depth", "nitrate", "phosphate", "silicate"])

nutrient_smoothed = interpolate_profiles(nutrient_binned, z_min = 3, z_max = 95, n_points = 400)

# 5) Shrink to 6m
small_smooth = nutrient_smoothed[nutrient_smoothed['depth']<95].copy()
small_smooth['depth'] = (small_smooth['depth'] / (95/6.21))-3.2

# -----------------------------
# Helpers for 3D advection
# -----------------------------

def _validate_columns(df: pd.DataFrame, required_cols: list, df_name: str):
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"{df_name} is missing required columns: {missing}")

def _interp_1d(x_src: np.ndarray, y_src: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    """Safely linearly interpolate y(x) to new x points. Extrapolation is clamped."""
    order = np.argsort(x_src)
    xs = np.asarray(x_src)[order]
    ys = np.asarray(y_src)[order]
    uniq_mask = np.r_[True, np.diff(xs) > 0]
    xs = xs[uniq_mask]
    ys = ys[uniq_mask]
    if xs.size == 1:
        return np.full_like(x_new, ys[0], dtype=float)
    x_new_clamped = np.clip(x_new, xs[0], xs[-1])
    return np.interp(x_new_clamped, xs, ys)

def _trilinear_interpolate(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    F: np.ndarray,
    xq: np.ndarray,
    yq: np.ndarray,
    zq: np.ndarray,
) -> np.ndarray:
    """Trilinear interpolation on a structured grid (x,y,z) for field F[z,y,x].
    All coordinate arrays must be strictly increasing. Queries are clamped inside bounds.
    Returns array with same shape as xq/yq/zq.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    # Clamp queries inside domain
    xq = np.clip(xq, x[0], x[-1])
    yq = np.clip(yq, y[0], y[-1])
    zq = np.clip(zq, z[0], z[-1])

    nx = x.size; ny = y.size; nz = z.size

    # Find left indices
    ix = np.searchsorted(x, xq, side='right') - 1
    iy = np.searchsorted(y, yq, side='right') - 1
    iz = np.searchsorted(z, zq, side='right') - 1

    ix = np.clip(ix, 0, nx - 2)
    iy = np.clip(iy, 0, ny - 2)
    iz = np.clip(iz, 0, nz - 2)

    # Compute normalized distances within the cell
    x0 = x[ix]; x1 = x[ix + 1]
    y0 = y[iy]; y1 = y[iy + 1]
    z0 = z[iz]; z1 = z[iz + 1]

    # Avoid divide-by-zero (degenerated spacing)
    dx = np.where(x1 > x0, x1 - x0, 1.0)
    dy = np.where(y1 > y0, y1 - y0, 1.0)
    dz = np.where(z1 > z0, z1 - z0, 1.0)

    tx = (xq - x0) / dx
    ty = (yq - y0) / dy
    tz = (zq - z0) / dz

    # Gather corner values (C-order: F[z,y,x])
    flat_shape = xq.shape
    ix_f = ix.ravel(); iy_f = iy.ravel(); iz_f = iz.ravel()

    def G(koff, joff, ioff):
        return F[ (iz_f + koff), (iy_f + joff), (ix_f + ioff) ]

    c000 = G(0, 0, 0)
    c001 = G(0, 0, 1)
    c010 = G(0, 1, 0)
    c011 = G(0, 1, 1)
    c100 = G(1, 0, 0)
    c101 = G(1, 0, 1)
    c110 = G(1, 1, 0)
    c111 = G(1, 1, 1)

    txf = tx.ravel(); tyf = ty.ravel(); tzf = tz.ravel()

    # Trilinear blend
    c00 = c000 * (1 - txf) + c001 * txf
    c01 = c010 * (1 - txf) + c011 * txf
    c10 = c100 * (1 - txf) + c101 * txf
    c11 = c110 * (1 - txf) + c111 * txf

    c0 = c00 * (1 - tyf) + c01 * tyf
    c1 = c10 * (1 - tyf) + c11 * tyf

    c = c0 * (1 - tzf) + c1 * tzf
    return c.reshape(flat_shape)

def _fill_missing_velocity_columns_2d(
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    nan_velocity_fill: str = "zero",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fill columns (vertical profiles) that are entirely NaN at a given (y,x) location
    by copying the nearest non-NaN column in the horizontal (x,y) plane.
    If nan_velocity_fill == 'nearest', also attempt nearest interpolation in-plane for any residual NaNs.
    Otherwise, set remaining NaNs to 0.

    u, v, w have shape (nz, ny, nx).
    """
    nz, ny, nx = w.shape

    # Identify columns fully NaN for all z
    col_nan = np.all(np.isnan(w), axis=0)  # shape (ny, nx)

    # If any columns are fully NaN, find nearest non-NaN donor
    if np.any(col_nan):
        Xg, Yg = np.meshgrid(x, y)
        coords = np.stack([Xg, Yg], axis=-1)  # (ny, nx, 2)
        valid = ~col_nan
        valid_idxs = np.argwhere(valid)
        if valid_idxs.size == 0:
            # If everything is NaN, zero-fill
            u[:] = np.nan_to_num(u, nan=0.0)
            v[:] = np.nan_to_num(v, nan=0.0)
            w[:] = np.nan_to_num(w, nan=0.0)
            return u, v, w

        donor_coords = coords[valid]  # (n_valid, 2)
        for (jy, ix) in np.argwhere(col_nan):
            # Find nearest donor column by Euclidean distance in (x,y)
            dxy = donor_coords - coords[jy, ix]
            dist2 = np.sum(dxy * dxy, axis=1)
            k = np.argmin(dist2)
            dyx = valid_idxs[k]  # (jy_d, ix_d)
            jy_d, ix_d = int(dyx[0]), int(dyx[1])
            u[:, jy, ix] = u[:, jy_d, ix_d]
            v[:, jy, ix] = v[:, jy_d, ix_d]
            w[:, jy, ix] = w[:, jy_d, ix_d]

    # Handle residual NaNs: vertical nearest interpolation within each column
    # (If scattered NaNs remain due to sparse z coverage.)
    for jy in range(ny):
        for ix in range(nx):
            if np.any(np.isnan(u[:, jy, ix])):
                col = u[:, jy, ix]
                good = ~np.isnan(col)
                if good.any():
                    idx = np.where(good, np.arange(len(col)), -1)
                    for k in range(1, len(col)):  # forward fill
                        if idx[k] == -1:
                            idx[k] = idx[k-1]
                    for k in range(len(col)-2, -1, -1):  # backward fill
                        if idx[k] == -1:
                            idx[k] = idx[k+1]
                    u[:, jy, ix] = col[idx]
            if np.any(np.isnan(v[:, jy, ix])):
                col = v[:, jy, ix]
                good = ~np.isnan(col)
                if good.any():
                    idx = np.where(good, np.arange(len(col)), -1)
                    for k in range(1, len(col)):
                        if idx[k] == -1:
                            idx[k] = idx[k-1]
                    for k in range(len(col)-2, -1, -1):
                        if idx[k] == -1:
                            idx[k] = idx[k+1]
                    v[:, jy, ix] = col[idx]
            if np.any(np.isnan(w[:, jy, ix])):
                col = w[:, jy, ix]
                good = ~np.isnan(col)
                if good.any():
                    idx = np.where(good, np.arange(len(col)), -1)
                    for k in range(1, len(col)):
                        if idx[k] == -1:
                            idx[k] = idx[k-1]
                    for k in range(len(col)-2, -1, -1):
                        if idx[k] == -1:
                            idx[k] = idx[k+1]
                    w[:, jy, ix] = col[idx]

    # Final NaN policy
    if nan_velocity_fill != 'nearest':
        u[:] = np.nan_to_num(u, nan=0.0)
        v[:] = np.nan_to_num(v, nan=0.0)
        w[:] = np.nan_to_num(w, nan=0.0)

    return u, v, w

def build_grid_and_initial_field(
    small_smooth: pd.DataFrame,
    velocity_example_file: str,
    nutrient_col: str = "nitrate",
    z_range=(-3.0, 3.0),
    dx=0.10,
    dy=0.10,
    dtype=np.float32,
    nan_velocity_fill="zero"
):
    """
    Build only:
      - x, y, z grids
      - initial nutrient 3D field
      - fast binning structure for velocity reuse

    No advection performed.
    """

    # -------------------------------------
    # Load sample velocity file to build grid
    # -------------------------------------
    df = pd.read_csv(velocity_example_file)

    required = [
        "Points:0", "Points:1", "Points:2",
    ]
    _validate_columns(df, required, "velocity example")

    x = pd.to_numeric(df["Points:0"], errors="coerce").to_numpy(dtype)
    y = pd.to_numeric(df["Points:1"], errors="coerce").to_numpy(dtype)
    z = pd.to_numeric(df["Points:2"], errors="coerce").to_numpy(dtype)

    base_mask = np.isfinite(x) & np.isfinite(y)
    x = x[base_mask]
    y = y[base_mask]
    z = z[base_mask]

    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())

    # -------------------------------------
    # Build X,Y regular grid via binning
    # -------------------------------------
    ix = np.rint((x - xmin) / dx).astype(np.int32)
    iy = np.rint((y - ymin) / dy).astype(np.int32)

    ix_bins = np.unique(ix)
    iy_bins = np.unique(iy)

    x_unique = (xmin + ix_bins * dx).astype(dtype)
    y_unique = (ymin + iy_bins * dy).astype(dtype)
    x_unique.sort()
    y_unique.sort()

    nx = len(x_unique)
    ny = len(y_unique)

    ix_to_axis = {b: i for i, b in enumerate(ix_bins)}
    iy_to_axis = {b: i for i, b in enumerate(iy_bins)}

    # -------------------------------------
    # Build vertical z grid from profile
    # -------------------------------------
    prof = small_smooth[["depth", nutrient_col]].dropna()
    prof["depth"] = pd.to_numeric(prof["depth"], errors="coerce")
    prof[nutrient_col] = pd.to_numeric(prof[nutrient_col], errors="coerce")

    prof = prof[(prof["depth"] >= z_range[0]) & (prof["depth"] <= z_range[1])]
    z_grid = np.sort(prof["depth"].to_numpy(dtype=float)).astype(dtype)

    nz = len(z_grid)

    # -------------------------------------
    # Build initial nutrient 3D field
    # -------------------------------------
    nutr0_1d = _interp_1d(
        prof["depth"].to_numpy(float),
        prof[nutrient_col].to_numpy(float),
        z_grid.astype(float)
    ).astype(dtype)

    nutrient_initial = nutr0_1d[:, None, None] * np.ones((1, ny, nx), dtype=dtype)

    # -------------------------------------
    # Precompute bin grouping (critical!)
    # -------------------------------------
    keys = (ix.astype(np.int64) << 32) ^ iy.astype(np.int64)
    order = np.argsort(keys, kind="mergesort")
    keys_sorted = keys[order]
    boundaries = np.flatnonzero(np.diff(keys_sorted) != 0) + 1
    starts = np.r_[0, boundaries]
    ends = np.r_[boundaries, len(keys_sorted)]

    return {
        "x": x_unique,
        "y": y_unique,
        "z": z_grid,
        "nutrient_initial": nutrient_initial,

        # Save for fast GIF
        "ix": ix,
        "iy": iy,
        "ix_to_axis": ix_to_axis,
        "iy_to_axis": iy_to_axis,
        "order": order,
        "starts": starts,
        "ends": ends,
        "base_mask": base_mask,
    }


def fast_load_velocity_file(
    vel_path: str,
    grid: dict,
    nan_velocity_fill="zero"
):
    """
    Fastest possible velocity loader:
    - Uses precomputed ix/iy groups
    - Loads CSV once
    - Interpolates into (nz, ny, nx)
    """

    df = pd.read_csv(vel_path)

    # Auto-detect velocity columns
    cols = df.columns
    possibilities = [
        ("U_average:0","U_average:1","U_average:2"),
        ("U:0","U:1","U:2"),
        ("velocity:0","velocity:1","velocity:2"),
        ("Velocity:0","Velocity:1","Velocity:2")
    ]
    for a,b,c in possibilities:
        if a in cols and b in cols and c in cols:
            ux, uy, uz = a, b, c
            break
    else:
        # fallback suffix match
        ux = [c for c in cols if c.endswith(":0")][0]
        uy = [c for c in cols if c.endswith(":1")][0]
        uz = [c for c in cols if c.endswith(":2")][0]

    # Unpack grid + binning
    x = grid["x"]
    y = grid["y"]
    z = grid["z"]
    ix = grid["ix"]
    iy = grid["iy"]
    order = grid["order"]
    starts = grid["starts"]
    ends = grid["ends"]
    ix_to_axis = grid["ix_to_axis"]
    iy_to_axis = grid["iy_to_axis"]
    base_mask = grid["base_mask"]

    nz=len(z); ny=len(y); nx=len(x)

    xv = df["Points:0"].to_numpy(float)[base_mask]
    yv = df["Points:1"].to_numpy(float)[base_mask]
    zv = df["Points:2"].to_numpy(float)[base_mask]

    u = df[ux].to_numpy(float)[base_mask]
    v = df[uy].to_numpy(float)[base_mask]
    w = df[uz].to_numpy(float)[base_mask]

    u3 = np.full((nz,ny,nx), np.nan)
    v3 = np.full((nz,ny,nx), np.nan)
    w3 = np.full((nz,ny,nx), np.nan)

    # groupwise interpolation
    for s,e in zip(starts, ends):

        idxs = order[s:e]
        ix_bin = ix[idxs[0]]
        iy_bin = iy[idxs[0]]

        ixa = ix_to_axis[ix_bin]
        iya = iy_to_axis[iy_bin]

        z_s = zv[idxs]
        u_s = u[idxs]
        v_s = v[idxs]
        w_s = w[idxs]

        good = np.isfinite(z_s)&np.isfinite(u_s)&np.isfinite(v_s)&np.isfinite(w_s)
        if not good.any():
            continue

        u3[:,iya,ixa] = _interp_1d(z_s[good],u_s[good],z)
        v3[:,iya,ixa] = _interp_1d(z_s[good],v_s[good],z)
        w3[:,iya,ixa] = _interp_1d(z_s[good],w_s[good],z)

    return _fill_missing_velocity_columns_2d(
        u3,v3,w3, x, y, nan_velocity_fill
    )


def animate_advection_y_slice(
    grid: dict,
    velocity_folder: str,
    y_value: float,
    total_time_s=150,
    dt_s=1.0,
    cmap="viridis",
    filename="y_slice.gif",
    fps=20,
    rectangle=None,
    progress=True,
    units_label="Nitrate (mmol m$^{-3}$)"
):
    """
    Clean, fast, self-contained advection animator.
    Uses fast_load_velocity_file() for each frame.
    """

    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.colors import Normalize
    from matplotlib.animation import FuncAnimation, PillowWriter
    import sys

    # Unpack grid + initial field
    x = grid["x"]; y = grid["y"]; z = grid["z"]
    C = grid["nutrient_initial"].copy()

    nz,ny,nx = C.shape
    iy_slice = np.argmin(np.abs(y - y_value))

    # Time stepping
    n_steps = int(np.ceil(total_time_s/dt_s))
    dt_eff = total_time_s/n_steps

    X, Z = np.meshgrid(x, z, indexing='xy')
    Xg,Yg = np.meshgrid(x,y, indexing='xy')

    X3 = np.tile(Xg,(nz,1,1))
    Y3 = np.tile(Yg,(nz,1,1))
    Z3 = np.tile(z.reshape(nz,1,1),(1,ny,nx))

    xmin,xmax = x[0],x[-1]
    ymin,ymax = y[0],y[-1]
    zmin,zmax = z[0],z[-1]

    # Setup plot
    vmin,vmax = float(C.min()), float(C.max())
    fig,ax = plt.subplots(figsize=(7,4),constrained_layout=True)
    norm = Normalize(vmin=vmin,vmax=vmax)

    im=ax.pcolormesh(X, Z, C[:,iy_slice,:], cmap=cmap, norm=norm, shading='auto')
    ax.set_xlabel("x-axis (m)")
    ax.set_ylabel("Depth (m)")

    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm,cmap=cmap),
        ax=ax,label=units_label,pad=0.005
    )

    if rectangle:
        x0,z0,wrect,hrect=rectangle
        rect_patch=patches.Rectangle((x0,z0),wrect,hrect,
            fill=True,facecolor="black",edgecolor="black",linewidth=0.1
        )
        ax.add_patch(rect_patch)
    else:
        rect_patch=None

    ax.autoscale(False)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(3,-3)

    time_text=ax.text(
        0.98,0.05,"t = 0 s",
        transform=ax.transAxes, ha="right",va="bottom",
        color="white",fontsize=10,
        bbox=dict(facecolor="black",alpha=0.5,edgecolor="none")
    )

    # Progress bar
    def print_progress(i,total,width=30):
        if not progress: return
        frac=i/total
        filled=int(width*frac)
        bar="█"*filled+"-"*(width-filled)
        pct=int(frac*100)
        sys.stdout.write(f"\rAnimating |{bar}| {pct:3d}%")
        sys.stdout.flush()
        if i==total:
            sys.stdout.write("\n"); sys.stdout.flush()

    frame={'i':0}

    # --- UPDATE FUNCTION ---
    def update(_):
        nonlocal C, im
        i=frame['i']

        # Load correct velocity field
        vel_path = os.path.join(
            velocity_folder, f"datasetNS5_{min(i,150):06d}.csv"
        )
        u3d,v3d,w3d = fast_load_velocity_file(vel_path, grid)

        X_dep = np.clip(X3 - u3d*dt_eff, xmin, xmax)
        Y_dep = np.clip(Y3 - v3d*dt_eff, ymin, ymax)
        Z_dep = np.clip(Z3 - w3d*dt_eff, zmin, zmax)

        C = _trilinear_interpolate(x,y,z, C, X_dep,Y_dep,Z_dep)

        im.remove()
        im=ax.pcolormesh(X, Z, C[:,iy_slice,:], cmap=cmap, norm=norm, shading='auto')

        if rect_patch:
            rect_patch.set_zorder(im.get_zorder()+1)

        time_text.set_text(f"t = {int(i*dt_eff)} s")

        print_progress(i,n_steps)

        frame['i']+=1
        return (im,)

    anim = FuncAnimation(fig,update,frames=n_steps,blit=False,interval=100)
    anim.save(filename, writer=PillowWriter(fps=fps))
    print_progress(n_steps,n_steps)
    plt.close(fig)

    return filename


# -------------------------
# 1. Build grid + initial field
# -------------------------
grid = build_grid_and_initial_field(
    small_smooth,
    velocity_example_file="/data/data1/kiyan/ns5TimeStepData/datasetNS5_000001.csv",
    nutrient_col="nitrate"
)

grid_p = build_grid_and_initial_field(
    small_smooth,
    velocity_example_file="/data/data1/kiyan/ns5TimeStepData/datasetNS5_000001.csv",
    nutrient_col="phosphate"
)

grid_si = build_grid_and_initial_field(
    small_smooth,
    velocity_example_file="/data/data1/kiyan/ns5TimeStepData/datasetNS5_000001.csv",
    nutrient_col="silicate"
)

# -------------------------
# 2. Generate GIF
# -------------------------
'''
gif = animate_advection_y_slice(
    grid=grid,
    velocity_folder="/data/data1/kiyan/ns5TimeStepData/",
    y_value=15.0,
    fps=5,
    filename="nitrate_y15_vertical.gif",
    rectangle=(9.0, -3.0, 2.0, 6.0)
)

print('Saved gif', gif)

gif = animate_advection_y_slice(
    grid=grid,
    velocity_folder="/data/data1/kiyan/ns5TimeStepData/",
    y_value=14.0,
    fps=5,
    filename="nitrate_y14_vertical.gif",
    rectangle=(0.0, 0.0, 0.0, 0.0)
)

print('Saved gif',gif)

gif = animate_advection_y_slice(
        grid=grid,
        velocity_folder="/data/data1/kiyan/ns5TimeStepData/",
        y_value=16.0,
        fps=5,
        filename="nitrate_y16_vertical.gif",
        rectangle=(0,0,0,0)
)
print("Saved gif",gif)
'''

gif = animate_advection_y_slice(
        grid=grid_p,
        velocity_folder="/data/data1/kiyan/ns5TimeStepData/",
        y_value=15.0,
        fps=5,
        filename="phosphate_y15_vertical.gif",
        rectangle=(9.0, -3.0, 2.0, 6.0)
)

print("Saved gif",gif)

gif = animate_advection_y_slice(
        grid=grid_si,
        velocity_folder="/data/data1/kiyan/ns5TimeStepData/",
        y_value=15.0,
        fps=5,
        filename="silicate_y15_vertical.gif",
        rectangle=(9.0, -3.0, 2.0, 6.0)
        )

print("Saved gif",gif)

