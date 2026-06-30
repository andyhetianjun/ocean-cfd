import numpy as np
import os, re

GAUGE_X = [1.0, 2.5, 4.0, 6.0, 7.5, 9.0]
STILL_WATER = 0.4

def get_time_dirs(case='.'):
    dirs = []
    for d in os.listdir(case):
        try:
            t = float(d)
            if t > 0 and os.path.isfile(os.path.join(case, d, 'alpha.water')):
                dirs.append((t, d))
        except ValueError:
            pass
    return sorted(dirs)

def read_scalar_field(filepath):
    with open(filepath) as f:
        content = f.read()
    match = re.search(r'internalField\s+nonuniform List<scalar>\s+\d+\s*\(([^)]+)\)',
                      content, re.DOTALL)
    if match:
        return np.array(match.group(1).split(), dtype=float)
    return None

print("Reading mesh cell centres ...")
os.system("postProcess -func writeCellCentres -time 0 > /dev/null 2>&1")

with open('0/C') as f:
    content = f.read()

match = re.search(r'internalField\s+nonuniform List<vector>\s+\d+\s*\((.+?)\)\s*;',
                  content, re.DOTALL)
raw = match.group(1).strip()
# strip parentheses from each (x y z) tuple
raw = raw.replace('(', '').replace(')', '')
coords = np.array(raw.split(), dtype=float).reshape(-1, 3)
cx, cy, cz = coords[:,0], coords[:,1], coords[:,2]
print(f"  {len(cx)} cells parsed")

print("Finding gauge columns ...")
gauge_cells = []
for gx in GAUGE_X:
    mask = (np.abs(cx - gx) < 0.15) & (np.abs(cy) < 0.3)
    idx = np.where(mask)[0]
    idx = idx[np.argsort(cz[idx])]
    gauge_cells.append((idx, cz[idx]))
    print(f"  x={gx}: {len(idx)} cells, z={cz[idx].min():.2f}..{cz[idx].max():.2f}")

print("\nReading time directories ...")
time_dirs = get_time_dirs()
print(f"  {len(time_dirs)} time steps found")

times, eta_all = [], [[] for _ in range(6)]

for t, tdir in time_dirs:
    alpha = read_scalar_field(os.path.join(tdir, 'alpha.water'))
    if alpha is None:
        for i in range(6): eta_all[i].append(STILL_WATER)
        times.append(t)
        continue
    for i, (idx, zc) in enumerate(gauge_cells):
        a_col = alpha[idx]
        dz = np.diff(zc, prepend=zc[0]-(zc[1]-zc[0]))
        eta_all[i].append(np.sum(a_col * dz))
    times.append(t)
    if len(times) % 20 == 0:
        print(f"  processed t={t:.2f}s")

times = np.array(times)
eta_all = [np.array(e) for e in eta_all]

np.savetxt('postProcessing/eta_gauges.dat',
           np.column_stack([times] + eta_all),
           header='time eta_G1 eta_G2 eta_G3 eta_G4 eta_G5 eta_G6',
           fmt='%.6f')
print("\nSaved: postProcessing/eta_gauges.dat")

def Hwave(eta, t):
    m = (t>=1.0)&(t<=11.0)
    return eta[m].max()-eta[m].min() if m.sum()>0 else np.nan

print("\n" + "="*50)
Hs = [Hwave(eta_all[i], times) for i in range(6)]
for i,h in enumerate(Hs):
    print(f"  G{i+1} x={GAUGE_X[i]:4.1f}m ({'up' if i<3 else 'down'}stream): H = {h:.4f} m")
Hi, Ht = np.nanmean(Hs[:3]), np.nanmean(Hs[3:])
print(f"\nH incident    : {Hi:.4f} m")
print(f"H transmitted : {Ht:.4f} m")
print(f"Kt            : {Ht/Hi:.3f}")
print("="*50)
