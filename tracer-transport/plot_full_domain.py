import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load('velocity_full_domain_05m/vel_595.npz')
g, c, x, y = d['grid_vel'], d['counts'], d['x'], d['y']
den = c.sum(axis=0)
u = np.where(den > 0, (g[0]*c).sum(axis=0)/np.maximum(den,1), np.nan)
v = np.where(den > 0, (g[1]*c).sum(axis=0)/np.maximum(den,1), np.nan)
print(f'empty columns: {(den==0).mean():.4f}')
print(f'upstream u (x<30): {np.nanmean(u[:, x<30]):.3f}  expect 0.496')

fig, ax = plt.subplots(figsize=(14, 3.2), constrained_layout=True)
m = np.nanpercentile(np.abs(v), 99)
im = ax.pcolormesh(x, y, v, shading='auto', cmap='RdBu_r', vmin=-m, vmax=m)
ax.set_aspect('equal')
ax.add_patch(plt.Rectangle((20,30), 200, 60, fill=False, ec='k', lw=1.5, ls='--'))
ax.add_patch(plt.Circle((46,60), 4, fc='0.3', ec='k', lw=0.8))
cax = ax.inset_axes([1.02, 0, 0.015, 1], transform=ax.transAxes)
fig.colorbar(im, cax=cax, label='depth-averaged v (m/s)')
ax.set_xlabel('x (m)'); ax.set_ylabel('y (m)')
ax.set_title('cross-stream velocity, depth-averaged, t = 595 s')
fig.savefig('full_domain_velocity.png', dpi=150)
print('saved full_domain_velocity.png')
