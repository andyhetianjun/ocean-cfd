#!/usr/bin/env python
"""
2D Advection Solver — one step up from 1D.

Solves pure advection in two dimensions (no diffusion):

    dc/dt + u*dc/dx + v*dc/dy = 0

This is the EXACT same idea as 1D, just applied in two directions and summed.
If you understood the 1D solver, this is mechanical:
  - 1D had one velocity u and one derivative dc/dx
  - 2D has velocities (u,v) and derivatives (dc/dx, dc/dy), added together

Same upwind rule (take the derivative from the direction the flow comes from),
same CFL timestep restriction (now accounting for BOTH directions).

The only genuinely new wrinkle: the CFL condition in 2D combines both
velocity components. A conservative form is:
    |u|*dt/dx + |v|*dt/dy  <=  1
"""

import numpy as np


def upwind_derivative(c, vel, dx, axis):
    """
    Compute the upwind spatial derivative of c along one axis.

    This is the SAME logic as 1D, factored into a reusable helper so we
    can call it once per direction (x, then y). In 3D we'll call it a
    third time (z) — the pattern doesn't change.

    c    : concentration field (2D array)
    vel  : velocity component along this axis (same shape as c)
    dx   : grid spacing along this axis
    axis : which array axis this direction corresponds to (0 or 1)
    """
    # backward difference uses the "previous" cell along this axis:
    #   (c[i] - c[i-1]) / dx   -> correct when vel > 0
    # forward difference uses the "next" cell:
    #   (c[i+1] - c[i]) / dx   -> correct when vel < 0
    back = (c - np.roll(c, 1, axis=axis)) / dx
    fwd  = (np.roll(c, -1, axis=axis) - c) / dx
    return np.where(vel >= 0, back, fwd)


def advect_2d(c0, u, v, dx, dy, dt, n_steps):
    """
    March 2D advection forward using first-order upwind.

    c0      : initial concentration, shape (ny, nx)
    u       : x-velocity, scalar or shape (ny, nx)  [m/s]
    v       : y-velocity, scalar or shape (ny, nx)  [m/s]
    dx, dy  : grid spacing in x and y [m]
    dt      : timestep [s]
    n_steps : number of steps

    Returns c_history: shape (n_steps+1, ny, nx)

    NOTE on array indexing: we store the field as c[y, x] (row=y, col=x).
      - axis 0 is y  -> uses dy and velocity v
      - axis 1 is x  -> uses dx and velocity u
    This (row=y, col=x) convention matches how images/plots are laid out.
    """
    c = c0.copy().astype(float)
    ny, nx = c.shape

    u = np.broadcast_to(np.asarray(u, dtype=float), (ny, nx)).copy()
    v = np.broadcast_to(np.asarray(v, dtype=float), (ny, nx)).copy()

    # CFL check combining both directions
    C = np.max(np.abs(u)) * dt / dx + np.max(np.abs(v)) * dt / dy
    if C > 1.0:
        raise ValueError(
            f"CFL violated: combined Courant C={C:.3f} > 1. Reduce dt."
        )

    c_history = np.empty((n_steps + 1, ny, nx))
    c_history[0] = c

    for step in range(n_steps):
        # x-direction: velocity u, axis 1, spacing dx
        dcdx = upwind_derivative(c, u, dx, axis=1)
        # y-direction: velocity v, axis 0, spacing dy
        dcdy = upwind_derivative(c, v, dy, axis=0)

        # update: dc/dt = -(u*dc/dx + v*dc/dy)
        c = c - dt * (u * dcdx + v * dcdy)

        c_history[step + 1] = c

    return c_history


# ---------------------------------------------------------------------------
# Test / demonstration
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # 2D domain: 10m x 10m, 150x150 cells
    Lx = Ly = 10.0
    nx = ny = 150
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # meshgrid so we can build a 2D initial blob. indexing='xy' gives shape (ny,nx).
    X, Y = np.meshgrid(x, y)   # both shape (ny, nx)

    # Uniform diagonal flow: 1 m/s in x, 0.5 m/s in y (blob should drift up-right)
    u = 1.0
    v = 0.5

    # CFL timestep for combined directions
    target_courant = 0.4
    dt = target_courant / (abs(u)/dx + abs(v)/dy)

    # Initial 2D Gaussian blob centered at (2, 2)
    c0 = np.exp(-(((X - 2.0) ** 2) + ((Y - 2.0) ** 2)) / (2 * 0.4 ** 2))

    # March long enough to drift ~5m in x (and ~2.5m in y since v=u/2)
    T = 5.0
    n_steps = int(T / dt)

    print(f"Domain: {Lx}x{Ly} m, grid {nx}x{ny}, dx={dx:.4f} dy={dy:.4f}")
    print(f"Velocity (u,v)=({u},{v}) m/s, dt={dt:.4f} s")
    print(f"Combined Courant={abs(u)*dt/dx + abs(v)*dt/dy:.3f}")
    print(f"Marching {n_steps} steps to T={T} s")
    print(f"Blob should drift: x by u*T={u*T:.1f}m, y by v*T={v*T:.1f}m")
    print(f"  from (2.0, 2.0) to ({2+u*T:.1f}, {2+v*T:.1f})")

    hist = advect_2d(c0, u, v, dx, dy, dt, n_steps)

    # find peak location at start and end (2D argmax)
    j0, i0 = np.unravel_index(np.argmax(hist[0]), hist[0].shape)
    j1, i1 = np.unravel_index(np.argmax(hist[-1]), hist[-1].shape)
    print(f"\nPeak start: (x={x[i0]:.2f}, y={y[j0]:.2f})")
    print(f"Peak end:   (x={x[i1]:.2f}, y={y[j1]:.2f})")
    print(f"Drift: x={x[i1]-x[i0]:.2f}m (expect {u*T:.1f}), "
          f"y={y[j1]-y[j0]:.2f}m (expect {v*T:.1f})")

    # mass conservation
    m0 = np.sum(hist[0]) * dx * dy
    m1 = np.sum(hist[-1]) * dx * dy
    print(f"\nTotal tracer: start={m0:.4f}, end={m1:.4f} "
          f"({100*(m1-m0)/m0:+.1f}% change)")
    print(f"Peak amplitude: {hist[0].max():.3f} -> {hist[-1].max():.3f} "
          f"(smearing expected from 1st-order upwind)")
