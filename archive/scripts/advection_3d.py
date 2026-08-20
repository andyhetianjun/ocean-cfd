#!/usr/bin/env python
"""
3D Advection Solver — the full advection equation.

Solves pure advection in three dimensions (no diffusion yet):

    dc/dt + u*dc/dx + v*dc/dy + w*dc/dz = 0

This IS the advection-only equation your manager gave you. It's the same
pattern as 1D and 2D, now in three directions:
  - 1D: one direction  (x)
  - 2D: two directions (x, y)
  - 3D: three directions (x, y, z)  <- here

The upwind helper is UNCHANGED from 2D — we just call it a third time for z.
That's the whole difference. Everything you validated in 1D/2D carries over.

CFL in 3D combines all three directions:
    |u|*dt/dx + |v|*dt/dy + |w|*dt/dz  <=  1
"""

import numpy as np


def upwind_derivative(c, vel, dx, axis):
    """
    Upwind spatial derivative of c along one axis. IDENTICAL to the 2D version.
    Called once per direction. In 3D: three calls (x, y, z).
    """
    back = (c - np.roll(c, 1, axis=axis)) / dx
    fwd  = (np.roll(c, -1, axis=axis) - c) / dx
    return np.where(vel >= 0, back, fwd)


def advect_3d(c0, u, v, w, dx, dy, dz, dt, n_steps, save_every=1):
    """
    March 3D advection forward using first-order upwind.

    c0         : initial concentration, shape (nz, ny, nx)
    u, v, w    : velocity components, scalar or shape (nz, ny, nx)  [m/s]
    dx, dy, dz : grid spacing in each direction [m]
    dt         : timestep [s]
    n_steps    : number of steps
    save_every : store a frame every N steps (3D histories are BIG in memory,
                 so we don't keep every single step by default)

    Returns
    -------
    c_history : shape (n_saved, nz, ny, nx) — saved frames
    saved_steps : the step index of each saved frame

    Array convention: c[z, y, x]
      axis 0 = z  -> dz, velocity w
      axis 1 = y  -> dy, velocity v
      axis 2 = x  -> dx, velocity u
    """
    c = c0.copy().astype(float)
    nz, ny, nx = c.shape

    u = np.broadcast_to(np.asarray(u, dtype=float), (nz, ny, nx)).copy()
    v = np.broadcast_to(np.asarray(v, dtype=float), (nz, ny, nx)).copy()
    w = np.broadcast_to(np.asarray(w, dtype=float), (nz, ny, nx)).copy()

    # CFL check combining all three directions
    C = (np.max(np.abs(u)) * dt / dx +
         np.max(np.abs(v)) * dt / dy +
         np.max(np.abs(w)) * dt / dz)
    if C > 1.0:
        raise ValueError(f"CFL violated: combined Courant C={C:.3f} > 1. Reduce dt.")

    frames = [c.copy()]
    saved_steps = [0]

    for step in range(1, n_steps + 1):
        dcdx = upwind_derivative(c, u, dx, axis=2)   # x
        dcdy = upwind_derivative(c, v, dy, axis=1)   # y
        dcdz = upwind_derivative(c, w, dz, axis=0)   # z

        # update: dc/dt = -(u*dcdx + v*dcdy + w*dcdz)
        c = c - dt * (u * dcdx + v * dcdy + w * dcdz)

        if step % save_every == 0 or step == n_steps:
            frames.append(c.copy())
            saved_steps.append(step)

    return np.array(frames), np.array(saved_steps)


# ---------------------------------------------------------------------------
# Test / demonstration
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # 3D domain: 10 x 10 x 10 m, 80^3 cells (kept modest so it runs fast)
    Lx = Ly = Lz = 10.0
    nx = ny = nz = 80
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    z = np.linspace(0, Lz, nz)
    dx = x[1]-x[0]; dy = y[1]-y[0]; dz = z[1]-z[0]

    # 3D meshgrid, indexing='ij' with (z,y,x) order gives arrays shaped (nz,ny,nx).
    Z, Y, X = np.meshgrid(z, y, x, indexing='ij')   # each shape (nz,ny,nx)

    # Uniform 3D flow: drift in all three directions
    u, v, w = 1.0, 0.5, 0.25

    target_courant = 0.4
    dt = target_courant / (abs(u)/dx + abs(v)/dy + abs(w)/dz)

    # Initial 3D Gaussian blob centered at (2,2,2)
    c0 = np.exp(-(((X-2)**2)+((Y-2)**2)+((Z-2)**2)) / (2*0.5**2))

    T = 4.0
    n_steps = int(T / dt)

    print(f"Domain {Lx}x{Ly}x{Lz} m, grid {nx}x{ny}x{nz} = {nx*ny*nz:,} cells")
    print(f"Velocity (u,v,w)=({u},{v},{w}) m/s, dt={dt:.4f} s")
    print(f"Combined Courant={abs(u)*dt/dx+abs(v)*dt/dy+abs(w)*dt/dz:.3f}")
    print(f"Marching {n_steps} steps to T={T} s")
    print(f"Blob should drift to ({2+u*T:.1f}, {2+v*T:.1f}, {2+w*T:.1f})")

    import time
    t0 = time.time()
    hist, steps = advect_3d(c0, u, v, w, dx, dy, dz, dt, n_steps, save_every=max(1,n_steps//10))
    elapsed = time.time() - t0

    # peak location start vs end
    k0,j0,i0 = np.unravel_index(np.argmax(hist[0]), hist[0].shape)
    k1,j1,i1 = np.unravel_index(np.argmax(hist[-1]), hist[-1].shape)
    print(f"\nPeak start: ({x[i0]:.2f}, {y[j0]:.2f}, {z[k0]:.2f})")
    print(f"Peak end:   ({x[i1]:.2f}, {y[j1]:.2f}, {z[k1]:.2f})")
    print(f"Drift: x={x[i1]-x[i0]:.2f} (exp {u*T:.1f}), "
          f"y={y[j1]-y[j0]:.2f} (exp {v*T:.1f}), "
          f"z={z[k1]-z[k0]:.2f} (exp {w*T:.1f})")

    m0 = np.sum(hist[0])*dx*dy*dz
    m1 = np.sum(hist[-1])*dx*dy*dz
    print(f"\nTotal tracer: start={m0:.4f}, end={m1:.4f} ({100*(m1-m0)/m0:+.1f}%)")
    print(f"Peak amplitude: {hist[0].max():.3f} -> {hist[-1].max():.3f}")
    print(f"\nCompute time: {elapsed:.2f} s for {n_steps} steps on {nx*ny*nz:,} cells")
