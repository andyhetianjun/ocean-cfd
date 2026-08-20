#!/usr/bin/env python
"""
Velocity Loader with CACHING + PARALLEL reading.

readmesh (computing cell centres) is the bottleneck (~20s/processor, mostly
parallelizable CPU work). Parallel across cores should cut this near-linearly
when the disk isn't contended. Read once, cache to .npz, use forever.
"""

import os
import glob
import time
import numpy as np
from fluidfoam import readvector, readmesh


def count_processors(case_dir):
    return len(glob.glob(os.path.join(case_dir, "processor*")))


def load_one_processor(args):
    case_dir, proc_id, timestep = args
    proc_dir = os.path.join(case_dir, f"processor{proc_id}")
    x, y, z = readmesh(proc_dir, verbose=False)
    U = readvector(proc_dir, timestep, "U", verbose=False)
    return np.column_stack([x, y, z]), U.T


def load_serial(case_dir, timestep, n_procs, verbose=True):
    all_pts, all_vel = [], []
    for p in range(n_procs):
        t0 = time.time()
        pts, vel = load_one_processor((case_dir, p, timestep))
        all_pts.append(pts); all_vel.append(vel)
        if verbose:
            print(f"  [{p+1}/{n_procs}] {pts.shape[0]} cells ({time.time()-t0:.1f}s)", flush=True)
    return np.vstack(all_pts), np.vstack(all_vel)


def load_parallel(case_dir, timestep, n_procs, workers, verbose=True):
    from multiprocessing import Pool
    args = [(case_dir, p, timestep) for p in range(n_procs)]
    all_pts, all_vel = [], []
    t0 = time.time()
    with Pool(workers) as pool:
        for i, (pts, vel) in enumerate(pool.imap(load_one_processor, args)):
            all_pts.append(pts); all_vel.append(vel)
            if verbose and (i % 16 == 0 or i == n_procs-1):
                print(f"  [{i+1}/{n_procs}] done ({time.time()-t0:.0f}s elapsed)", flush=True)
    return np.vstack(all_pts), np.vstack(all_vel)


def main(case_dir, timestep, out_path, workers=1, max_procs=None):
    n = count_processors(case_dir)
    if max_procs:
        n = min(n, max_procs)
    print(f"Processors: {n}, timestep: {timestep}, workers: {workers}\n", flush=True)

    t0 = time.time()
    if workers > 1:
        points, velocity = load_parallel(case_dir, timestep, n, workers)
    else:
        points, velocity = load_serial(case_dir, timestep, n)
    read_time = time.time() - t0

    print(f"\nRead {points.shape[0]:,} cells in {read_time/60:.2f} min ({read_time:.0f}s)", flush=True)
    print(f"Domain: x[{points[:,0].min():.1f},{points[:,0].max():.1f}] "
          f"y[{points[:,1].min():.1f},{points[:,1].max():.1f}] "
          f"z[{points[:,2].min():.1f},{points[:,2].max():.1f}]")
    Umag = np.sqrt((velocity**2).sum(axis=1))
    print(f"|U|: min={Umag.min():.3f} max={Umag.max():.3f} mean={Umag.mean():.3f}")

    # only cache when reading the FULL set (not a partial test)
    if max_procs is None:
        np.savez_compressed(out_path, points=points, velocity=velocity, timestep=timestep)
        print(f"\nCACHED to {out_path}.npz ({os.path.getsize(out_path+'.npz')/1e6:.0f} MB)")
    else:
        print(f"\n(partial test of {n} procs — not caching)")


if __name__ == "__main__":
    import sys
    base = "/shared_folder/yongxing/OpenFOAM/simulationCases/flow_past_cylinder/realistic_cases/domain_D8X600Y120/uniform"
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    max_procs = int(sys.argv[2]) if len(sys.argv) > 2 else None
    main(base, "450", "velocity_450_cache", workers=workers, max_procs=max_procs)
