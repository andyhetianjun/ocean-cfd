#!/usr/bin/env python
"""
Velocity Loader for decomposed OpenFOAM cases.

Reads a velocity field from an OpenFOAM case decomposed across many
processor folders (processor0, processor1, ...), WITHOUT reconstructPar.
Each processor holds a chunk of the domain; we read them all and merge into
a single point cloud, then (later) bin onto a regular grid for the solver.

Reusable: works on any decomposed case, any timestep. Nothing hardcoded.
"""

import os
import glob
import numpy as np
from fluidfoam import readvector, readmesh


def count_processors(case_dir):
    return len(glob.glob(os.path.join(case_dir, "processor*")))


def load_one_processor(args):
    """Read cell-center coords + velocity for one processor. Top-level
    function so it works with multiprocessing later."""
    case_dir, proc_id, timestep = args
    proc_dir = os.path.join(case_dir, f"processor{proc_id}")
    x, y, z = readmesh(proc_dir, verbose=False)
    U = readvector(proc_dir, timestep, "U", verbose=False)
    pts = np.column_stack([x, y, z])
    vel = U.T
    return pts, vel


def load_all_processors(case_dir, timestep, n_procs=None, verbose=True):
    """Read ALL processors (serial) and merge into one point cloud."""
    if n_procs is None:
        n_procs = count_processors(case_dir)
    all_pts, all_vel = [], []
    for p in range(n_procs):
        pts, vel = load_one_processor((case_dir, p, timestep))
        all_pts.append(pts)
        all_vel.append(vel)
        if verbose and (p % 32 == 0 or p == n_procs - 1):
            print(f"  read processor {p}/{n_procs-1} ({pts.shape[0]} cells)")
    return np.vstack(all_pts), np.vstack(all_vel)


if __name__ == "__main__":
    import time
    base = "/shared_folder/yongxing/OpenFOAM/simulationCases/flow_past_cylinder/realistic_cases/domain_D8X600Y120/uniform"

    n = count_processors(base)
    print(f"Case has {n} processors")
    print(f"Reading ALL {n} processors at timestep 450 (serial)...\n")

    t0 = time.time()
    points, velocity = load_all_processors(base, "450")
    elapsed = time.time() - t0

    print(f"\n=== MERGED FULL DOMAIN ===")
    print(f"Total cells: {points.shape[0]:,}")
    print(f"Read time: {elapsed:.1f} s (serial)")
    print(f"\nFull domain bounds:")
    print(f"  x: {points[:,0].min():.2f} to {points[:,0].max():.2f}")
    print(f"  y: {points[:,1].min():.2f} to {points[:,1].max():.2f}")
    print(f"  z: {points[:,2].min():.2f} to {points[:,2].max():.2f}")

    Umag = np.sqrt((velocity**2).sum(axis=1))
    print(f"\nVelocity magnitude across FULL domain:")
    print(f"  min={Umag.min():.3f}, max={Umag.max():.3f}, mean={Umag.mean():.3f}")
    print(f"  (max >> mean would indicate wake/vortex structure is present)")

    imax = np.argmax(Umag)
    print(f"\nFastest flow at: ({points[imax,0]:.1f}, {points[imax,1]:.1f}, "
          f"{points[imax,2]:.1f}) with |U|={Umag[imax]:.3f}")
