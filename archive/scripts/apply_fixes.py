import ast, shutil, sys
SRC = "oil_spill.py"
src = open(SRC).read()

def sub(old, new, what):
    global src
    n = src.count(old)
    if n != 1: sys.exit("ABORT: %s found %d times, expected 1" % (what, n))
    src = src.replace(old, new)

# 1. mass total per progress print (needs re-advection to appear)
sub('            print(f"    step {n+1}/{len(TIMESTEPS)}", flush=True)',
'''            _dv = float(x[1]-x[0]) * float(y[1]-y[0]) * float(z_tr[1]-z_tr[0])
            print(f"    step {n+1}/{len(TIMESTEPS)}  mass={C.sum()*_dv/1000.0:.4f} kg",
                  flush=True)''', 'mass diagnostic')

# 2. units: g/m3 -> g/m2 at plot time, case 1 only
sub('    lab = "Oil (g m$^{-2}$)"; nut = "Oil"; fp = f"oil_{BEARING}"',
'''    lab = "Oil (g m$^{-2}$)" if not VERTICAL_MIX else "Oil (g m$^{-3}$)"
    nut = "Oil"; fp = f"oil_{BEARING}"''', 'label')

sub('    for (mode, sv, fn, ttl), slabs in zip(jobs, all_slabs):',
'''    if not VERTICAL_MIX:
        # case 1: all oil in one layer, so g/m3 * dz is exact areal density.
        # case 2 spreads vertically and needs a column integral, not a slice.
        _dz = float(z_tr[1] - z_tr[0])
        all_slabs = [[s * _dz for s in sl] for sl in all_slabs]
    for (mode, sv, fn, ttl), slabs in zip(jobs, all_slabs):''', 'unit scaling')

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_fixes")
open(SRC, "w").write(src)
print("OK")
