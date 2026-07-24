import numpy as np
d = np.load('velocity_full_domain_05m/vel_595.npz')
c, x, y = d['counts'], d['x'], d['y']
f = c.sum(axis=0) > 0
px = f[(y >= 40) & (y < 80), :].mean(axis=0)
py = f[:, x < 250].mean(axis=1)
print('x profile (within y 40-80):')
for lo in range(0, 600, 40):
    m = (x >= lo) & (x < lo + 40)
    print('  x %3d-%3d: %.3f' % (lo, lo + 40, px[m].mean()))
print('y profile (within x < 250):')
for lo in range(0, 120, 10):
    m = (y >= lo) & (y < lo + 10)
    print('  y %3d-%3d: %.3f' % (lo, lo + 10, py[m].mean()))
for p, ax, n in ((px, x, 'x'), (py, y, 'y')):
    h = p > 0.9
    if h.any():
        print('%s: coverage above 0.9 spans %.1f to %.1f' % (n, ax[h][0], ax[h][-1]))
