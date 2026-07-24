import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load('velocity_full_domain_05m/vel_595.npz')
c, x, y = d['counts'], d['x'], d['y']
f = (c.sum(axis=0) > 0).astype(np.float32)
B = 8
ny, nx = f.shape
fb = f.reshape(ny//B, B, nx//B, B).mean(axis=(1,3))
xe = np.linspace(x[0]-0.25, x[-1]+0.25, nx//B + 1)
ye = np.linspace(y[0]-0.25, y[-1]+0.25, ny//B + 1)

fig, ax = plt.subplots(figsize=(14, 3.2), constrained_layout=True)
im = ax.pcolormesh(xe, ye, fb, cmap='viridis', vmin=0, vmax=1)
ax.set_aspect('equal')
ax.add_patch(plt.Rectangle((20,30), 200, 60, fill=False, ec='w', lw=1.6, ls='--'))
ax.add_patch(plt.Circle((46,60), 4, fc='w', ec='k', lw=0.8))
cax = ax.inset_axes([1.02, 0, 0.015, 1], transform=ax.transAxes)
fig.colorbar(im, cax=cax, label='fraction of 0.5 m bins containing mesh points')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_title('LES mesh coverage over full 600 x 120 m domain (dashed = tracer window)')
fig.savefig('mesh_coverage.png', dpi=150)
print('saved mesh_coverage.png')
