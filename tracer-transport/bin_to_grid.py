#!/usr/bin/env python
"""
Bin a scattered velocity point cloud onto a regular grid.

LES gives velocity at ~44M scattered cell centers. The finite-difference
solver needs a REGULAR grid. This bins scattered data onto that grid by
averaging points per cell. Reusable: any point cloud, region, resolution.
"""

import numpy as np


def bin_to_grid(points, velocity, bounds, shape):
    """
    Bin scattered velocity onto a regular grid by averaging.
    points   : (N,3) coords;  velocity : (N,3) (u,v,w)
    bounds   : ((xmin,xmax),(ymin,ymax),(zmin,zmax))
    shape    : (nx, ny, nz)
    Returns grid_vel (3,nz,ny,nx) [NaN where empty], (x,y,z) coords, counts.
    """
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounds
    nx, ny, nz = shape

    m = ((points[:,0] >= xmin) & (points[:,0] < xmax) &
         (points[:,1] >= ymin) & (points[:,1] < ymax) &
         (points[:,2] >= zmin) & (points[:,2] < zmax))
    pts = points[m]; vel = velocity[m]

    ix = np.clip(((pts[:,0]-xmin)/(xmax-xmin)*nx).astype(int), 0, nx-1)
    iy = np.clip(((pts[:,1]-ymin)/(ymax-ymin)*ny).astype(int), 0, ny-1)
    iz = np.clip(((pts[:,2]-zmin)/(zmax-zmin)*nz).astype(int), 0, nz-1)
    flat = (iz * ny + iy) * nx + ix

    ncells = nx * ny * nz
    counts = np.bincount(flat, minlength=ncells).astype(float)
    grid_vel = np.full((3, ncells), np.nan)
    for c in range(3):
        s = np.bincount(flat, weights=vel[:,c], minlength=ncells)
        with np.errstate(invalid='ignore', divide='ignore'):
            grid_vel[c] = np.where(counts > 0, s / counts, np.nan)

    grid_vel = grid_vel.reshape(3, nz, ny, nx)
    counts = counts.reshape(nz, ny, nx)
    x = xmin + (np.arange(nx)+0.5)*(xmax-xmin)/nx
    y = ymin + (np.arange(ny)+0.5)*(ymax-ymin)/ny
    z = zmin + (np.arange(nz)+0.5)*(zmax-zmin)/nz
    return grid_vel, (x, y, z), counts


if __name__ == "__main__":
    import time
    print("Loading cached velocity...", flush=True)
    d = np.load("velocity_450_cache.npz")
    points, velocity = d["points"], d["velocity"]
    print(f"  {points.shape[0]:,} scattered cells loaded", flush=True)

    # Focus region: cylinder (~46m) + downstream wake, full width & depth.
    # x from 20 (just upstream) to 220 (near+mid wake); y,z full domain.
    bounds = ((20.0, 220.0), (0.0, 119.4), (-19.8, 19.8))
    # Resolution: ~1m cells in x, ~2m in y, ~2m in z -> 200 x 60 x 20
    shape = (200, 60, 20)   # nx, ny, nz

    print(f"\nBinning to grid {shape} over region {bounds}...", flush=True)
    t0 = time.time()
    gvel, (x, y, z), counts = bin_to_grid(points, velocity, bounds, shape)
    print(f"  done in {time.time()-t0:.1f}s", flush=True)

    filled = (counts > 0).sum()
    print(f"\nGrid cells with data: {filled:,}/{counts.size:,} "
          f"({100*filled/counts.size:.0f}%)")
    print(f"Empty cells (cylinder interior / sparse): {counts.size-filled:,}")

    gmag = np.sqrt(np.nansum(gvel**2, axis=0))
    print(f"Gridded |U|: min={np.nanmin(gmag):.3f} "
          f"max={np.nanmax(gmag):.3f} mean={np.nanmean(gmag):.3f}")

    # save the gridded velocity for the solver
    np.savez_compressed("velocity_450_gridded",
                        grid_vel=gvel, x=x, y=y, z=z, counts=counts,
                        bounds=np.array(bounds), shape=np.array(shape))
    print(f"\nSaved velocity_450_gridded.npz — ready for the solver")
