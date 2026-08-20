#!/usr/bin/env python
"""Tracer gifs -- fixed color scaling so tracer is vivid."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation, PillowWriter

CYL_X, CYL_Y, CYL_R = 46.0, 60.0, 4.0


def make_vertical(frames, x, y, z, saved_t, out="tracer_vertical.gif", fps=8):
    jy = np.argmin(np.abs(y - CYL_Y))
    slabs = frames[:, :, jy, :]
    # fixed vmax at a fraction of the peak so tracer is vivid, not washed out
    vmax = float(np.percentile(frames[frames>0.01], 95)) if (frames>0.01).any() else 1.0
    norm = Normalize(0, vmax)
    X, Zg = np.meshgrid(x, z)
    fig, ax = plt.subplots(figsize=(9,4), constrained_layout=True)
    im = ax.pcolormesh(X, Zg, slabs[0], cmap="viridis", norm=norm, shading="auto")
    ax.set_xlabel("x (m)"); ax.set_ylabel("Depth (m)"); ax.invert_yaxis()
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax,
                 label="Tracer concentration", pad=0.01)
    ax.add_patch(patches.Rectangle((CYL_X-CYL_R, z.min()), 2*CYL_R, z.max()-z.min(),
                                   facecolor="black", zorder=5))
    tt = ax.text(0.98,0.05,"",transform=ax.transAxes,ha="right",va="bottom",
                 color="white",fontsize=10,bbox=dict(facecolor="black",alpha=0.5,edgecolor="none"))
    def update(f):
        im.set_array(slabs[f].ravel()); tt.set_text(f"t = {saved_t[f]:.0f} s"); return im,tt
    FuncAnimation(fig,update,frames=len(slabs),blit=False).save(out,writer=PillowWriter(fps=fps))
    plt.close(fig); return out


def make_plan(frames, x, y, z, saved_t, out="tracer_planview.gif", fps=8):
    plan = frames.max(axis=1)
    vmax = float(np.percentile(frames[frames>0.01], 95)) if (frames>0.01).any() else 1.0
    norm = Normalize(0, vmax)
    X, Yg = np.meshgrid(x, y)
    fig, ax = plt.subplots(figsize=(10,4), constrained_layout=True)
    im = ax.pcolormesh(X, Yg, plan[0], cmap="viridis", norm=norm, shading="auto")
    ax.set_xlabel("x (m)  [downstream ->]"); ax.set_ylabel("y (m)  [width]")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax,
                 label="Tracer (max over depth)", pad=0.01)
    ax.add_patch(patches.Circle((CYL_X, CYL_Y), CYL_R, facecolor="black",
                                edgecolor="white", zorder=5))
    tt = ax.text(0.98,0.05,"",transform=ax.transAxes,ha="right",va="bottom",
                 color="white",fontsize=10,bbox=dict(facecolor="black",alpha=0.5,edgecolor="none"))
    def update(f):
        im.set_array(plan[f].ravel()); tt.set_text(f"t = {saved_t[f]:.0f} s"); return im,tt
    FuncAnimation(fig,update,frames=len(plan),blit=False).save(out,writer=PillowWriter(fps=fps))
    plt.close(fig); return out


if __name__ == "__main__":
    d = np.load("tracer_frames.npz")
    frames=d["frames"]; x,y,z=d["x"],d["y"],d["z"]; saved_t=d["saved_t"]
    print(f"Loaded {len(frames)} frames {frames.shape[1:]}", flush=True)
    print(f"  {make_vertical(frames,x,y,z,saved_t)}", flush=True)
    print(f"  {make_plan(frames,x,y,z,saved_t)}", flush=True)
    print("Done.", flush=True)
