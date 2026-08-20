#!/usr/bin/env python
"""Vertical diffusivity passes 1 and 2.  Usage: diffusivity.py {1|2}"""
import os, sys
import numpy as np

VEL_DIR   = "data/velocity_full_025"
TIMESTEPS = [str(t) for t in range(1049, 1401)]
DT_S      = 1.0
TRACER_NZ = 100

MODE_N  = 1        # cosine mode number
AMP_A   = 5.0      # amplitude
MEAN_CM = 10.0     # offset, must exceed AMP_A so C stays positive
TAG     = f"n{MODE_N}"

def load_velocity(t):
    d = np.load(os.path.join(VEL_DIR, f"vel_{t}.npz"))
    g = d["grid_vel"]
    return (np.nan_to_num(g[0]).astype(np.float32),
            np.nan_to_num(g[1]).astype(np.float32),
            np.nan_to_num(g[2]).astype(np.float32))

def load_grid():
    d = np.load(os.path.join(VEL_DIR, f"vel_{TIMESTEPS[0]}.npz"))
    x = d["x"].astype(float); y = d["y"].astype(float)
    z_vel = d["z"].astype(float)
    z_tr = np.linspace(z_vel[0], z_vel[-1], TRACER_NZ) if TRACER_NZ != len(z_vel) else z_vel
    return x, y, z_vel, z_tr

def velocity_on_tracer_z(u, v, w, z_vel, z_tr):
    if z_tr.shape == z_vel.shape and np.allclose(z_tr, z_vel):
        return u, v, w
    iz = np.clip(np.searchsorted(z_vel, z_tr, side="right") - 1, 0, len(z_vel) - 2)
    z0 = z_vel[iz]; z1 = z_vel[iz + 1]
    t = ((z_tr - z0) / np.where(z1 > z0, z1 - z0, 1.0))[:, None, None].astype(np.float32)
    lift = lambda F: F[iz] * (1 - t) + F[iz + 1] * t
    return lift(u), lift(v), lift(w)

def _trilinear_interpolate(x, y, z, F, xq, yq, zq):
    xq = np.clip(xq, x[0], x[-1]); yq = np.clip(yq, y[0], y[-1]); zq = np.clip(zq, z[0], z[-1])
    nx, ny, nz = x.size, y.size, z.size
    ix = np.clip(np.searchsorted(x, xq, side="right") - 1, 0, nx - 2)
    iy = np.clip(np.searchsorted(y, yq, side="right") - 1, 0, ny - 2)
    iz = np.clip(np.searchsorted(z, zq, side="right") - 1, 0, nz - 2)
    x0, x1 = x[ix], x[ix + 1]; y0, y1 = y[iy], y[iy + 1]; z0, z1 = z[iz], z[iz + 1]
    tx = ((xq - x0) / np.where(x1 > x0, x1 - x0, 1.0)).astype(np.float32)
    ty = ((yq - y0) / np.where(y1 > y0, y1 - y0, 1.0)).astype(np.float32)
    tz = ((zq - z0) / np.where(z1 > z0, z1 - z0, 1.0)).astype(np.float32)
    shape = xq.shape
    ixf, iyf, izf = ix.ravel(), iy.ravel(), iz.ravel()
    G = lambda k, j, i: F[izf + k, iyf + j, ixf + i]
    txf, tyf, tzf = tx.ravel(), ty.ravel(), tz.ravel()
    c00 = G(0,0,0)*(1-txf) + G(0,0,1)*txf
    c01 = G(0,1,0)*(1-txf) + G(0,1,1)*txf
    c10 = G(1,0,0)*(1-txf) + G(1,0,1)*txf
    c11 = G(1,1,0)*(1-txf) + G(1,1,1)*txf
    c0 = c00*(1-tyf) + c01*tyf
    c1 = c10*(1-tyf) + c11*tyf
    return (c0*(1-tzf) + c1*tzf).reshape(shape)

def cosine_field(z_tr, ny, nx):
    zb = float(z_tr[0]); H = float(z_tr[-1] - z_tr[0])
    prof = MEAN_CM + AMP_A * np.cos(MODE_N * np.pi * (z_tr - zb) / H)
    return np.repeat(np.repeat(prof.astype(np.float32)[:,None,None], ny, 1), nx, 2)

def run(which):
    x, y, z_vel, z_tr = load_grid()
    nz, ny, nx = len(z_tr), len(y), len(x)
    print(f"grid {nx} x {ny} x {nz}   mode n={MODE_N}  A={AMP_A}  Cm={MEAN_CM}", flush=True)
    C = cosine_field(z_tr, ny, nx)
    Z3, Y3, X3 = (a.astype(np.float32) for a in np.meshgrid(z_tr, y, x, indexing="ij"))
    xmin, xmax = x[0], x[-1]; ymin, ymax = y[0], y[-1]; zmin, zmax = z_tr[0], z_tr[-1]
    nt = len(TIMESTEPS)

    if which == 1:
        sum_c = np.zeros((nz, ny, nx), np.float64)
        sum_w = np.zeros((nz, ny, nx), np.float64)
    else:
        m = np.load(f"data/means_{TAG}.npz")
        cbar = m["cbar"]; wbar = m["wbar"]
        F_time  = np.zeros((nz, nt)); F_layer = np.zeros((nz, nt))
        Cl      = np.zeros((nz, nt)); mass    = np.zeros(nt)

    for n, t in enumerate(TIMESTEPS):
        u, v, w = load_velocity(t)
        u, v, w = velocity_on_tracer_z(u, v, w, z_vel, z_tr)
        if which == 1:
            sum_c += C; sum_w += w
        else:
            cp = C - cbar; wp = w - wbar
            F_time[:, n]  = (wp * cp).mean(axis=(1, 2))
            wl = w.mean(axis=(1, 2)); cl = C.mean(axis=(1, 2))
            F_layer[:, n] = (w * C).mean(axis=(1, 2)) - wl * cl
            Cl[:, n] = cl
            mass[n]  = float(C.sum())
        Xd = np.clip(X3 - u*np.float32(DT_S), xmin, xmax)
        Yd = np.clip(Y3 - v*np.float32(DT_S), ymin, ymax)
        Zd = np.clip(Z3 - w*np.float32(DT_S), zmin, zmax)
        C = _trilinear_interpolate(x, y, z_tr, C, Xd, Yd, Zd)
        del u, v, w, Xd, Yd, Zd
        if (n+1) % 20 == 0:
            print(f"    pass{which} step {n+1}/{nt}", flush=True)

    if which == 1:
        np.savez(f"data/means_{TAG}.npz",
                 cbar=(sum_c/nt).astype(np.float32),
                 wbar=(sum_w/nt).astype(np.float32), z=z_tr)
        print("wrote means", flush=True)
    else:
        np.savez(f"data/flux_{TAG}.npz", F_time=F_time, F_layer=F_layer,
                 Cl=Cl, mass=mass, z=z_tr, nt=nt)
        print("wrote flux", flush=True)

if __name__ == "__main__":
    run(int(sys.argv[1]))
