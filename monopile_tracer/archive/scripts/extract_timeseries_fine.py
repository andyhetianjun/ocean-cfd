#!/usr/bin/env python
"""
Multi-timestep velocity extraction -> small gridded cache. RESUMABLE.

The mesh is static (fixed cylinder), so cell centres are identical at every
timestep. We reuse the points already cached from t=450 and only read the U
field per timestep (0.4s vs 19.5s for readmesh). Bin indices precomputed once.
Each timestep saves a small gridded .npz. Re-running skips finished timesteps.
"""
import os, re, time
import numpy as np
from fluidfoam import readvector
from multiprocessing import Pool

CASE = "/shared_folder/yongxing/OpenFOAM/simulationCases/flow_past_cylinder/realistic_cases/domain_D8X600Y120/uniform"
POINTS_CACHE = "velocity_450_cache.npz"
OUT_DIR = "velocity_timeseries_fine"
BOUNDS = ((20.0, 220.0), (30.0, 90.0), (-20.0, 20.0))
SHAPE = (800, 240, 84)
WORKERS = 24
N_PROCS = 256
MAX_TIMESTEPS = 200


def discover_timesteps(case):
    p0 = os.path.join(case, "processor0")
    times = []
    for d in os.listdir(p0):
        full = os.path.join(p0, d)
        if os.path.isdir(full) and re.fullmatch(r"[0-9]+(\.[0-9]+)?", d):
            if os.path.exists(os.path.join(full, "U")):
                times.append(d)
    times.sort(key=float)
    return times


def read_U_one_proc(args):
    case, proc_id, t = args
    U = readvector(os.path.join(case, f"processor{proc_id}"), t, "U", verbose=False)
    return U.T


def precompute_bin_index(points, bounds, shape):
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounds
    nx, ny, nz = shape
    m = ((points[:,0] >= xmin) & (points[:,0] < xmax) &
         (points[:,1] >= ymin) & (points[:,1] < ymax) &
         (points[:,2] >= zmin) & (points[:,2] < zmax))
    p = points[m]
    ix = np.clip(((p[:,0]-xmin)/(xmax-xmin)*nx).astype(np.int32), 0, nx-1)
    iy = np.clip(((p[:,1]-ymin)/(ymax-ymin)*ny).astype(np.int32), 0, ny-1)
    iz = np.clip(((p[:,2]-zmin)/(zmax-zmin)*nz).astype(np.int32), 0, nz-1)
    flat = (iz.astype(np.int64)*ny + iy)*nx + ix
    counts = np.bincount(flat, minlength=nx*ny*nz).astype(np.float32)
    x = xmin + (np.arange(nx)+0.5)*(xmax-xmin)/nx
    y = ymin + (np.arange(ny)+0.5)*(ymax-ymin)/ny
    z = zmin + (np.arange(nz)+0.5)*(zmax-zmin)/nz
    return m, flat, counts, (x, y, z)


def bin_velocity(vel_in_region, flat, counts, shape):
    nx, ny, nz = shape
    ncells = nx*ny*nz
    out = np.full((3, ncells), np.nan, dtype=np.float32)
    for c in range(3):
        s = np.bincount(flat, weights=vel_in_region[:,c], minlength=ncells)
        with np.errstate(invalid='ignore', divide='ignore'):
            out[c] = np.where(counts > 0, s/np.maximum(counts,1e-12), np.nan)
    return out.reshape(3, nz, ny, nx)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Loading cached cell-centre points (mesh static, reused for all timesteps)...", flush=True)
    t0 = time.time()
    points = np.load(POINTS_CACHE)["points"]
    print(f"  {points.shape[0]:,} points in {time.time()-t0:.0f}s", flush=True)

    print("Precomputing bin indices once...", flush=True)
    t0 = time.time()
    mask, flat, counts, (gx, gy, gz) = precompute_bin_index(points, BOUNDS, SHAPE)
    print(f"  done in {time.time()-t0:.0f}s ({mask.sum():,} points in region)", flush=True)
    del points

    times = discover_timesteps(CASE)
    print(f"\nTimesteps with U available: {len(times)}", flush=True)
    print(f"  first: {times[:5]} ... last: {times[-5:]}", flush=True)

    todo = [t for t in times if float(t) > 0 and not os.path.exists(os.path.join(OUT_DIR, f"vel_{t}.npz"))]
    todo = todo[:MAX_TIMESTEPS]
    if not todo:
        print("Nothing to do."); return
    print(f"Doing {len(todo)} timesteps this run (resumable): {todo[0]} .. {todo[-1]}\n", flush=True)

    t_start = time.time()
    with Pool(WORKERS) as pool:
        for n, t in enumerate(todo, 1):
            ts = time.time()
            args = [(CASE, p, t) for p in range(N_PROCS)]
            chunks = pool.map(read_U_one_proc, args)
            vel = np.vstack(chunks); del chunks
            if vel.shape[0] != mask.shape[0]:
                print(f"  [skip] t={t}: got {vel.shape[0]} values, expected {mask.shape[0]}", flush=True)
                del vel; continue
            gvel = bin_velocity(vel[mask], flat, counts, SHAPE); del vel
            np.savez_compressed(os.path.join(OUT_DIR, f"vel_{t}"),
                                grid_vel=gvel, x=gx, y=gy, z=gz,
                                counts=counts.reshape(SHAPE[2],SHAPE[1],SHAPE[0]), timestep=t)
            gm = np.sqrt(np.nansum(gvel**2, axis=0))
            print(f"  [{n}/{len(todo)}] t={t}  {time.time()-ts:.0f}s  "
                  f"|U| max={np.nanmax(gm):.3f} mean={np.nanmean(gm):.3f}  "
                  f"(total {(time.time()-t_start)/60:.1f} min)", flush=True)

    print(f"\nDONE. {len(todo)} timesteps in {(time.time()-t_start)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
