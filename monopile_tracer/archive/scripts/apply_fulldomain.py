import ast, shutil, sys
SRC = "nutrient_tracer_fmt.py"
src = open(SRC).read()
lines = src.split("\n")

def find(sub):
    h = [i for i,l in enumerate(lines) if sub in l]
    if len(h) != 1: sys.exit("ABORT: %d hits for %r" % (len(h), sub))
    return h[0]

# 1. config flag near top -- insert after the TRACER_NZ / VEL_DIR block by
#    placing it right before build_initial_field def
i = find("def build_initial_field")
cfg = [
 "FULL_DOMAIN = True          # pad computed 20-220 window out to full 0-600 for plotting",
 "FULL_X = (0.0, 600.0)       # full LES domain extent in x",
 "",
]
lines[i:i] = cfg

# 2. padding helper: inserted just before def main
i = find("def main():")
helper = '''def pad_to_full(slabs, x, C0, mode, slice_idx):
    """Pad each computed frame in x out to FULL_X, filling outside the
    computed window with the t=0 initial profile value (per depth row for
    xz, single plane value for xy). Returns (padded_slabs, x_full)."""
    dx = x[1] - x[0]
    xf0 = FULL_X[0] + dx/2.0
    nx_full = int(round((FULL_X[1] - FULL_X[0]) / dx))
    x_full = xf0 + np.arange(nx_full) * dx
    j0 = int(round((x[0] - x_full[0]) / dx))          # where computed window starts
    if mode == "xz":
        fill_col = C0[:, 0, 0][:, None]               # (nz,1) initial value per depth
    else:
        fill_col = C0[slice_idx, :, 0][:, None]       # (ny,1) plane value per y-row
    out = []
    for s in slabs:
        base = np.repeat(fill_col, nx_full, axis=1).astype(np.float32)
        base[:, j0:j0 + s.shape[1]] = s
        out.append(base)
    return out, x_full

'''
lines[i:i] = helper.split("\n")

src = "\n".join(lines)

# 3. widen plan-view figsize for the 5:1 shape (xz keeps (8,4))
src = src.replace(
    'figsize, xl, yl = (9, 4), "x (m)", "y (m)"',
    'figsize, xl, yl = ((14, 3.5) if FULL_DOMAIN else (9, 4)), "x (m)", "y (m)"')

# 4. in main: build z-plane index map and pad each slab before animate
src = src.replace(
    '    for (mode, sv, fn, ttl), slabs in zip(jobs, all_slabs):\n'
    '        t0 = time.time()\n'
    '        print(f"  {fn}", flush=True)\n'
    '        animate(mode, sv, slabs, x, y, z_tr, fn, title_str=ttl, units_label=lab)',
    '    for (mode, sv, fn, ttl), slabs in zip(jobs, all_slabs):\n'
    '        t0 = time.time()\n'
    '        print(f"  {fn}", flush=True)\n'
    '        xplot = x\n'
    '        if FULL_DOMAIN:\n'
    '            siz = int(np.argmin(np.abs(z_tr - sv))) if mode == "xy" else 0\n'
    '            slabs, xplot = pad_to_full(slabs, x, C0, mode, siz)\n'
    '        animate(mode, sv, slabs, xplot, y, z_tr, fn, title_str=ttl, units_label=lab)')

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_fd")
open(SRC, "w").write(src)
print("OK, backup nutrient_tracer_fmt.py.bak_fd")
