#!/usr/bin/env python
"""
1D Advection Solver — the foundation.

Solves the pure advection equation (no diffusion):

    dc/dt + u * dc/dx = 0

which says: a concentration c is carried along by a velocity u without
spreading. A blob of tracer just drifts; its shape is preserved (ideally).

This is phase 1 of the tracer tool: advection only. We get the numerics
right here in 1D — where we can SEE the answer — before scaling to 3D.

Two numerical ideas that MUST be right, or it explodes:

1. UPWIND DIFFERENCING for the space derivative dc/dx.
   The naive "centered" difference (c[i+1]-c[i-1])/(2dx) causes oscillations
   that overshoot, go negative (unphysical for concentration), and blow up.
   The fix: take the derivative in the direction the flow is COMING FROM
   ("upwind"). If u>0 (flow moving right), the info comes from the left,
   so we use the backward difference (c[i]-c[i-1])/dx. If u<0, forward diff.

2. CFL CONDITION for the timestep.
   The flow must not cross more than one grid cell per timestep, or the
   scheme can't "see" where the tracer came from. Formally:
       Courant number  C = |u| * dt / dx  <=  1
   We pick dt to keep C at a safe value (e.g. 0.4). Violate this -> blows up.
"""

import numpy as np


def advect_1d(c0, u, dx, dt, n_steps):
    """
    March the 1D advection equation forward in time using first-order upwind.

    Parameters
    ----------
    c0      : initial concentration, shape (nx,)
    u       : velocity, either a scalar or an array shape (nx,)  [m/s]
    dx      : grid spacing [m]
    dt      : timestep [s]
    n_steps : how many timesteps to march

    Returns
    -------
    c_history : shape (n_steps+1, nx) — concentration at every timestep
                (we keep all frames so we can animate/plot the evolution)
    """
    c = c0.copy().astype(float)
    nx = c.size

    # allow u to be a scalar (uniform flow) or a per-cell array (varying flow)
    u = np.broadcast_to(np.asarray(u, dtype=float), (nx,)).copy()

    # sanity check the CFL condition up front, so we fail loudly not silently
    C = np.max(np.abs(u)) * dt / dx
    if C > 1.0:
        raise ValueError(
            f"CFL violated: Courant number C={C:.3f} > 1. "
            f"Reduce dt below {dx/np.max(np.abs(u)):.4g} s, or coarsen the flow."
        )

    c_history = np.empty((n_steps + 1, nx))
    c_history[0] = c

    for step in range(n_steps):
        # We compute the spatial derivative dc/dx with upwind differencing.
        #
        # backward difference: (c[i] - c[i-1]) / dx   -> used where u > 0
        # forward  difference: (c[i+1] - c[i]) / dx   -> used where u < 0
        #
        # np.roll shifts the array so we can do this vectorized (no python loop
        # over cells). roll(c, 1) puts c[i-1] at position i; roll(c,-1) puts c[i+1].
        dcdx_backward = (c - np.roll(c, 1)) / dx    # uses left neighbor
        dcdx_forward  = (np.roll(c, -1) - c) / dx   # uses right neighbor

        # pick the upwind derivative cell-by-cell based on the sign of u:
        #   where u >= 0 use backward, where u < 0 use forward
        dcdx = np.where(u >= 0, dcdx_backward, dcdx_forward)

        # the update: dc/dt = -u * dc/dx  ->  c_new = c - dt * u * dc/dx
        c = c - dt * u * dcdx

        c_history[step + 1] = c

    return c_history


# ---------------------------------------------------------------------------
# Test / demonstration
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Build a 1D domain: 0 to 10 meters, 200 cells
    L = 10.0
    nx = 200
    x = np.linspace(0, L, nx)
    dx = x[1] - x[0]

    # Uniform flow to the right at 1 m/s
    u = 1.0

    # Pick timestep from the CFL condition: aim for Courant ~ 0.4 (safe margin)
    target_courant = 0.4
    dt = target_courant * dx / abs(u)

    # Initial condition: a Gaussian "blob" of tracer centered at x=2
    c0 = np.exp(-((x - 2.0) ** 2) / (2 * 0.3 ** 2))

    # March forward long enough for the blob to travel a few meters.
    # In time T, the blob moves u*T meters. We want it to move ~5 m, so T=5 s.
    T = 5.0
    n_steps = int(T / dt)

    print(f"Domain: 0..{L} m, nx={nx}, dx={dx:.4f} m")
    print(f"Velocity u={u} m/s, dt={dt:.4f} s, Courant={abs(u)*dt/dx:.3f}")
    print(f"Marching {n_steps} steps to T={T} s")
    print(f"Blob should move u*T = {u*T:.1f} m: from x=2.0 to x={2.0+u*T:.1f}")

    hist = advect_1d(c0, u, dx, dt, n_steps)

    # Check: where is the peak at start vs end?
    peak_start = x[np.argmax(hist[0])]
    peak_end   = x[np.argmax(hist[-1])]
    print(f"\nPeak started at x={peak_start:.2f}, ended at x={peak_end:.2f}")
    print(f"Moved {peak_end - peak_start:.2f} m (expected ~{u*T:.1f} m)")

    # Check mass conservation (total tracer should be ~preserved)
    mass_start = np.sum(hist[0]) * dx
    mass_end   = np.sum(hist[-1]) * dx
    print(f"\nTotal tracer: start={mass_start:.4f}, end={mass_end:.4f} "
          f"({100*(mass_end-mass_start)/mass_start:+.1f}% change)")

    # Check peak amplitude (upwind diffuses/smears the peak — this is expected)
    print(f"Peak amplitude: start={hist[0].max():.4f}, end={hist[-1].max():.4f} "
          f"(smearing from numerical diffusion is normal for 1st-order upwind)")
