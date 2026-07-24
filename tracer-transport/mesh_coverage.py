import numpy as np
d = np.load('velocity_full_domain_05m/vel_595.npz')
c, x, y = d['counts'], d['x'], d['y']
den = c.sum(axis=0); filled = den > 0
nzp = (c > 0).sum(axis=0)[filled]
print('z-bins populated per filled column: median %.0f of 84, min %d' % (np.median(nzp), nzp.min()))
print('overall filled: %.3f' % filled.mean())
for lo, hi in ((0,100),(100,200),(200,300),(300,400),(400,500),(500,600)):
    m = (x >= lo) & (x < hi)
    print('  x %3d-%3d: filled %.3f' % (lo, hi, filled[:, m].mean()))
for lo, hi in ((0,30),(30,50),(50,70),(70,90),(90,120)):
    m = (y >= lo) & (y < hi)
    print('  y %3d-%3d: filled %.3f' % (lo, hi, filled[m, :].mean()))
