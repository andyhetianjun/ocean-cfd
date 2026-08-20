#!/usr/bin/env python
"""Vorticity-Z animation, three panels.  Run from the folder containing data/."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches, gridspec
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation, PillowWriter

VEL_DIR   = "data/velocity_full_025"
TIMESTEPS = [str(t) for t in range(1049, 1401)]
DT_S      = 1.0
FPS       = 24
OUTFILE   = "vorticity_z_3panel.gif"
CYL_X, CYL_Y, CYL_R = 40.0, 60.0, 4.0
PLAN_Z    = 0.0
SECT_X    = 100.0
SECT_Y    = 60.0
X_MAX     = 300.0
VLIM      = 0.15
DX        = 0.25

def load_slices(t, iz, ix, iy):
    d = np.load(os.path.join(VEL_DIR, f"vel_{t}.npz"))
    g = d["grid_vel"]
    u = np.nan_to_num(g[0]).astype(np.float32)
    v = np.nan_to_num(g[1]).astype(np.float32)
    del g, d
    plan = np.gradient(v[iz], DX, axis=1) - np.gradient(u[iz], DX, axis=0)
    cross = (v[:, :, ix+1] - v[:, :, ix-1])/(2*DX) - np.gradient(u[:, :, ix], DX, axis=1)
    stream = np.gradient(v[:, iy, :], DX, axis=1) - (u[:, iy+1, :] - u[:, iy-1, :])/(2*DX)
    del u, v
    return plan, cross, stream

def main():
    d0 = np.load(os.path.join(VEL_DIR, f"vel_{TIMESTEPS[0]}.npz"))
    x = d0["x"].astype(float); y = d0["y"].astype(float); z = d0["z"].astype(float)
    del d0
    iz = int(np.argmin(np.abs(z - PLAN_Z)))
    ix = int(np.argmin(np.abs(x - SECT_X)))
    iy = int(np.argmin(np.abs(y - SECT_Y)))
    nx_crop = int(np.searchsorted(x, X_MAX))
    print(f"grid {len(x)} x {len(y)} x {len(z)}", flush=True)
    print(f"plan z={z[iz]:.2f} (i={iz})  cross x={x[ix]:.2f} (i={ix})  stream y={y[iy]:.2f} (i={iy})", flush=True)
    print(f"crop x to {x[nx_crop-1]:.1f} m, vlim +/-{VLIM}", flush=True)

    xc = x[:nx_crop]
    norm = Normalize(vmin=-VLIM, vmax=VLIM); cmap = "RdBu_r"
    fig = plt.figure(figsize=(15, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[3,1], height_ratios=[1,1],
                           hspace=0.28, wspace=0.18, left=0.06, right=0.90,
                           top=0.93, bottom=0.07)
    ax_plan = fig.add_subplot(gs[0,0]); ax_cross = fig.add_subplot(gs[0,1])
    ax_stream = fig.add_subplot(gs[1,0])
    plan, cross, stream = load_slices(TIMESTEPS[0], iz, ix, iy)

    im_plan = ax_plan.pcolormesh(xc, y, plan[:, :nx_crop], cmap=cmap, norm=norm, shading="auto")
    ax_plan.add_patch(patches.Circle((CYL_X, CYL_Y), CYL_R, facecolor="white",
                                     edgecolor="black", lw=0.8, zorder=5))
    ax_plan.axvline(SECT_X, color="0.3", lw=0.7); ax_plan.axhline(SECT_Y, color="0.3", lw=0.7)
    ax_plan.set_xlabel("x (m)"); ax_plan.set_ylabel("y (m)")
    ax_plan.set_title(f"z = {z[iz]:.1f} m (mid-depth)", loc="left", fontsize=10)
    ax_plan.set_aspect(1.0)

    im_cross = ax_cross.pcolormesh(z, y, cross.T, cmap=cmap, norm=norm, shading="auto")
    ax_cross.set_xlabel("z (m)"); ax_cross.set_ylabel("y (m)")
    ax_cross.set_title(f"x = {SECT_X:.0f} m", loc="left", fontsize=10)

    im_stream = ax_stream.pcolormesh(xc, z, stream[:, :nx_crop], cmap=cmap, norm=norm, shading="auto")
    ax_stream.axvline(CYL_X, color="0.3", lw=0.7)
    ax_stream.set_xlabel("x (m)"); ax_stream.set_ylabel("z (m)")
    ax_stream.set_title(f"y = {SECT_Y:.0f} m (pile centreline)", loc="left", fontsize=10)

    cax = fig.add_axes([0.92, 0.30, 0.015, 0.40])
    fig.colorbar(im_plan, cax=cax, label="Vorticity Z (1/s)", extend="both")
    tt = fig.text(0.50, 0.965, "", ha="center", fontsize=13, fontweight="bold")
    fig.text(0.92, 0.20, "Cylinder: D = 8 m\nInlet: Ux = 0.5 m/s\n"
             "Domain: 600 x 120 x 40 m\nPile at (40, 60)",
             fontsize=8, va="top", family="monospace")

    def update(i):
        p, c, s = load_slices(TIMESTEPS[i], iz, ix, iy)
        im_plan.set_array(p[:, :nx_crop].ravel())
        im_cross.set_array(c.T.ravel())
        im_stream.set_array(s[:, :nx_crop].ravel())
        tt.set_text(f"Vorticity Z    Time: {int(TIMESTEPS[0]) + int(i*DT_S)} s")
        if (i+1) % 20 == 0:
            print(f"    frame {i+1}/{len(TIMESTEPS)}", flush=True)
        return (im_plan, im_cross, im_stream, tt)

    print("rendering...", flush=True)
    FuncAnimation(fig, update, frames=len(TIMESTEPS), blit=False).save(
        os.path.join("results", OUTFILE), writer=PillowWriter(fps=FPS))
    plt.close(fig)
    print(f"saved results/{OUTFILE}", flush=True)

if __name__ == "__main__":
    main()
