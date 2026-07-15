#!/usr/bin/env python
"""Run tracer through LES wake -- with cylinder mask + continuous source."""
import numpy as np


def upwind_derivative(c, vel, dx, axis):
    back = (c - np.roll(c, 1, axis=axis)) / dx
    fwd  = (np.roll(c, -1, axis=axis) - c) / dx
    return np.where(vel >= 0, back, fwd)


def advect_step(c, u, v, w, dx, dy, dz, dt):
    dcdx = upwind_derivative(c, u, dx, axis=2)
    dcdy = upwind_derivative(c, v, dy, axis=1)
    dcdz = upwind_derivative(c, w, dz, axis=0)
    return c - dt * (u*dcdx + v*dcdy + w*dcdz)


def main():
    d = np.load("velocity_450_gridded.npz")
    gvel = d["grid_vel"]; x,y,z = d["x"],d["y"],d["z"]; counts = d["counts"]
    u,v,w = np.nan_to_num(gvel[0]), np.nan_to_num(gvel[1]), np.nan_to_num(gvel[2])
    nz,ny,nx = u.shape
    dx=x[1]-x[0]; dy=y[1]-y[0]; dz=z[1]-z[0]

    solid = (counts == 0)                # cylinder = cells with no fluid
    print(f"Solid (cylinder) cells: {solid.sum()}", flush=True)

    umax=max(abs(u).max(),1e-9); vmax=max(abs(v).max(),1e-9); wmax=max(abs(w).max(),1e-9)
    dt = 0.4 / (umax/dx + vmax/dy + wmax/dz)

    # continuous source upstream of cylinder, centered in y,z
    Z,Y,X = np.meshgrid(z,y,x, indexing='ij')
    source = np.exp(-(((X-25)**2)/(2*2**2)+((Y-60)**2)/(2*6**2)+((Z-0)**2)/(2*6**2)))
    source[source < 0.05] = 0.0

    c = np.zeros_like(u)
    T = 400.0
    n_steps = int(T/dt)
    print(f"grid {nx}x{ny}x{nz}, dt={dt:.3f}, steps={n_steps}, T={T}s", flush=True)
    print(f"max vel u={umax:.2f} v={vmax:.2f} w={wmax:.2f}", flush=True)

    import time; t0=time.time()
    frames=[]; saved_t=[]
    for step in range(n_steps):
        c = np.maximum(c, source)       # steady release
        c = advect_step(c,u,v,w,dx,dy,dz,dt)
        c[solid] = 0.0                  # enforce solid cylinder
        c = np.clip(c, 0, None)
        if step % max(1,n_steps//40) == 0:
            frames.append(c.copy()); saved_t.append(step*dt)
    print(f"advected in {time.time()-t0:.1f}s", flush=True)

    frames=np.array(frames)
    prof = frames[-1].max(axis=(0,1))
    reached = x[prof > 0.05]
    print(f"peak conc: {frames[-1].max():.3f}", flush=True)
    if len(reached): print(f"tracer reached x={reached.max():.0f} (source x=25)", flush=True)
    np.savez_compressed("tracer_frames", frames=frames, x=x,y=y,z=z, saved_t=np.array(saved_t))
    print(f"saved {len(frames)} frames", flush=True)


if __name__ == "__main__":
    main()
