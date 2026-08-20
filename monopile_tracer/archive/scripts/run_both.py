#!/usr/bin/env python
"""Run BOTH pulse and continuous releases through the LES wake, make gifs."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation, PillowWriter

CYL_X, CYL_Y, CYL_R = 46.0, 60.0, 4.0


def upwind_derivative(c, vel, dx, axis):
    back = (c - np.roll(c, 1, axis=axis)) / dx
    fwd  = (np.roll(c, -1, axis=axis) - c) / dx
    return np.where(vel >= 0, back, fwd)


def advect_step(c, u, v, w, dx, dy, dz, dt):
    dcdx = upwind_derivative(c, u, dx, axis=2)
    dcdy = upwind_derivative(c, v, dy, axis=1)
    dcdz = upwind_derivative(c, w, dz, axis=0)
    return c - dt * (u*dcdx + v*dcdy + w*dcdz)


def run(mode, u, v, w, x, y, z, solid, dx, dy, dz, dt, T):
    Z, Y, X = np.meshgrid(z, y, x, indexing='ij')
    blob = np.exp(-(((X-25)**2)/(2*2**2)+((Y-60)**2)/(2*6**2)+((Z-0)**2)/(2*6**2)))
    blob[blob < 0.05] = 0.0
    n_steps = int(T/dt)
    c = blob.copy() if mode == "pulse" else np.zeros_like(u)
    frames, saved_t = [], []
    for step in range(n_steps):
        if mode == "continuous":
            c = np.maximum(c, blob)
        c = advect_step(c, u, v, w, dx, dy, dz, dt)
        c[solid] = 0.0
        c = np.clip(c, 0, None)
        if step % max(1, n_steps//40) == 0:
            frames.append(c.copy()); saved_t.append(step*dt)
    return np.array(frames), np.array(saved_t)


def make_gif(frames, x, y, z, saved_t, view, out, fps=8):
    vmax = float(np.percentile(frames[frames>0.01], 95)) if (frames>0.01).any() else 1.0
    norm = Normalize(0, vmax)
    if view == "vertical":
        jy = np.argmin(np.abs(y - CYL_Y)); data = frames[:, :, jy, :]
        H, V = np.meshgrid(x, z); figsize, xl, yl, invert = (9,4), "x (m)", "Depth (m)", True
    else:
        data = frames.max(axis=1)
        H, V = np.meshgrid(x, y); figsize, xl, yl, invert = (10,4), "x (m)  [downstream ->]", "y (m)  [width]", False
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.pcolormesh(H, V, data[0], cmap="viridis", norm=norm, shading="auto")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    if invert: ax.invert_yaxis()
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax, label="Tracer concentration", pad=0.01)
    if view == "vertical":
        ax.add_patch(patches.Rectangle((CYL_X-CYL_R, z.min()), 2*CYL_R, z.max()-z.min(), facecolor="black", zorder=5))
    else:
        ax.add_patch(patches.Circle((CYL_X, CYL_Y), CYL_R, facecolor="black", edgecolor="white", zorder=5))
    tt = ax.text(0.98,0.05,"",transform=ax.transAxes,ha="right",va="bottom",color="white",fontsize=10,
                 bbox=dict(facecolor="black",alpha=0.5,edgecolor="none"))
    def upd(f):
        im.set_array(data[f].ravel()); tt.set_text(f"t = {saved_t[f]:.0f} s"); return im,tt
    FuncAnimation(fig,upd,frames=len(data),blit=False).save(out,writer=PillowWriter(fps=fps))
    plt.close(fig); return out


def main():
    d=np.load("velocity_450_gridded.npz")
    g=d["grid_vel"]; x,y,z=d["x"],d["y"],d["z"]; counts=d["counts"]
    u,v,w=np.nan_to_num(g[0]),np.nan_to_num(g[1]),np.nan_to_num(g[2])
    solid=(counts==0)
    nz,ny,nx=u.shape
    dx=x[1]-x[0]; dy=y[1]-y[0]; dz=z[1]-z[0]
    um=max(abs(u).max(),1e-9);vm=max(abs(v).max(),1e-9);wm=max(abs(w).max(),1e-9)
    dt=0.4/(um/dx+vm/dy+wm/dz)
    for mode in ["pulse","continuous"]:
        print(f"--- {mode} ---", flush=True)
        fr,st = run(mode,u,v,w,x,y,z,solid,dx,dy,dz,dt,T=400.0)
        for view in ["vertical","plan"]:
            out=f"tracer_{mode}_{view}.gif"
            make_gif(fr,x,y,z,st,view,out); print(f"  saved {out}", flush=True)


if __name__ == "__main__":
    main()
