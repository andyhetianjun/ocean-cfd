#!/usr/bin/env python
"""Vorticity-Z animation, two panels. Run from the folder containing data/."""
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
OUTFILE   = "speed_2panel.gif"
CYL_X, CYL_Y, CYL_R = 40.0, 60.0, 4.0
PLAN_Z    = 0.0
SECT_Y    = 60.0
X_MAX     = 300.0
VMIN      = 0.40      # m/s
VMAX      = 0.60      # m/s (freestream 0.5)
DX        = 0.25
Z_SURFACE = 19.761904761904763

def load_slices(t, iz, iy):
    d = np.load(os.path.join(VEL_DIR, f"vel_{t}.npz"))
    g3 = d["grid_vel"]
    u = np.nan_to_num(g3[0]).astype(np.float32)
    v = np.nan_to_num(g3[1]).astype(np.float32)
    w = np.nan_to_num(g3[2]).astype(np.float32)
    plan = np.sqrt(u[iz]**2 + v[iz]**2 + w[iz]**2)
    stream = np.sqrt(u[:, iy, :]**2 + v[:, iy, :]**2 + w[:, iy, :]**2)
    del u, v, w
    return plan, stream

def main():
    d0 = np.load(os.path.join(VEL_DIR, f"vel_{TIMESTEPS[0]}.npz"))
    x = d0["x"].astype(float); y = d0["y"].astype(float); z = d0["z"].astype(float)
    del d0
    iz = int(np.argmin(np.abs(z - PLAN_Z)))
    iy = int(np.argmin(np.abs(y - SECT_Y)))
    nx_crop = int(np.searchsorted(x, X_MAX))
    xc = x[:nx_crop]
    depth = z - Z_SURFACE

    print(f"grid {len(x)} x {len(y)} x {len(z)}", flush=True)
    print(f"plan z={z[iz]:.2f} (i={iz})   stream y={y[iy]:.2f} (i={iy})", flush=True)
    print(f"depth axis {depth[-1]:.1f} (surface) to {depth[0]:.1f} (bed)", flush=True)
    print(f"crop x to {xc[-1]:.1f} m, speed {VMIN}-{VMAX} m/s", flush=True)

    norm = Normalize(vmin=VMIN, vmax=VMAX); cmap = "viridis"
    fig = plt.figure(figsize=(13, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.0, 1.35], hspace=0.11,
                           left=0.08, right=0.86, top=0.92, bottom=0.09)
    ax_plan = fig.add_subplot(gs[0])
    ax_stream = fig.add_subplot(gs[1], sharex=ax_plan)

    plan, stream = load_slices(TIMESTEPS[0], iz, iy)

    im_plan = ax_plan.pcolormesh(xc, y, plan[:, :nx_crop], cmap=cmap, norm=norm, shading="auto")
    ax_plan.add_patch(patches.Circle((CYL_X, CYL_Y), CYL_R, facecolor="black",
                                     edgecolor="black", lw=0.8, zorder=5))
    ax_plan.axhline(SECT_Y, color="0.4", lw=0.6)
    ax_plan.set_ylabel("y (m)")
    ax_plan.set_title(f"z = {z[iz] - Z_SURFACE:.0f} m (mid-depth)", loc="left", fontsize=10)
    ax_plan.tick_params(labelbottom=False)

    im_stream = ax_stream.pcolormesh(xc, depth, stream[:, :nx_crop], cmap=cmap,
                                     norm=norm, shading="auto")
    ax_stream.add_patch(patches.Rectangle(
        (CYL_X - CYL_R, depth[0]), 2*CYL_R, depth[-1] - depth[0],
        facecolor="black", edgecolor="black", zorder=5))
    ax_stream.set_xlabel("x (m)"); ax_stream.set_ylabel("z (m)")
    ax_stream.set_title(f"y = {SECT_Y:.0f} m (pile centreline)", loc="left", fontsize=10)
    ax_stream.set_ylim(depth[0], depth[-1])
    ax_stream.set_yticks([0, -5, -10, -15, -20, -25, -30, -35, -40])
    ax_stream.set_xlim(xc[0], xc[-1])

    cax = fig.add_axes([0.88, 0.09, 0.018, 0.83])
    fig.colorbar(im_plan, cax=cax, label="Speed (m s$^{-1}$)", extend="both")
    tt = fig.text(0.47, 0.965, "", ha="center", fontsize=13, fontweight="bold")

    def update(i):
        p, s = load_slices(TIMESTEPS[i], iz, iy)
        im_plan.set_array(p[:, :nx_crop].ravel())
        im_stream.set_array(s[:, :nx_crop].ravel())
        tt.set_text(f"Speed    Time: {int(TIMESTEPS[0]) + int(i*DT_S)} s")
        if (i+1) % 20 == 0:
            print(f"    frame {i+1}/{len(TIMESTEPS)}", flush=True)
        return (im_plan, im_stream, tt)

    print("rendering...", flush=True)
    FuncAnimation(fig, update, frames=len(TIMESTEPS), blit=False).save(
        os.path.join("results", OUTFILE), writer=PillowWriter(fps=FPS))
    plt.close(fig)
    print(f"saved results/{OUTFILE}", flush=True)

if __name__ == "__main__":
    main()
