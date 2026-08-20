import ast, shutil, sys
SRC = "nutrient_tracer_fmt.py"
src = open(SRC).read()

# 1. add FULL_Y next to FULL_X
old_fx = 'FULL_X = (0.0, 600.0)       # full LES domain extent in x'
new_fx = old_fx + '\nFULL_Y = (0.0, 120.0)       # full LES domain extent in y'
assert src.count(old_fx) == 1, "FULL_X line not unique"
src = src.replace(old_fx, new_fx)

# 2. replace pad_to_full with x+y version
old_helper = '''def pad_to_full(slabs, x, C0, mode, slice_idx):
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
    return out, x_full'''

new_helper = '''def pad_to_full(slabs, x, y, C0, mode, slice_idx):
    """Pad each computed frame out to the full LES domain, filling outside
    the computed window with the t=0 initial profile value. x is padded for
    all modes; y is padded only for xy plan views (where the vertical axis
    is y). Returns (padded_slabs, x_full, y_full)."""
    dx = x[1] - x[0]
    xf0 = FULL_X[0] + dx/2.0
    nx_full = int(round((FULL_X[1] - FULL_X[0]) / dx))
    x_full = xf0 + np.arange(nx_full) * dx
    j0 = int(round((x[0] - x_full[0]) / dx))          # computed window start in x
    if mode == "xz":
        # rows are depth; fill each row with its initial value; no y padding
        fill_col = C0[:, 0, 0][:, None]               # (nz,1)
        out = []
        for s in slabs:
            base = np.repeat(fill_col, nx_full, axis=1).astype(np.float32)
            base[:, j0:j0 + s.shape[1]] = s
            out.append(base)
        return out, x_full, y
    # xy plan view: single plane value, pad both x and y
    plane_val = float(C0[slice_idx, 0, 0])
    dy = y[1] - y[0]
    yf0 = FULL_Y[0] + dy/2.0
    ny_full = int(round((FULL_Y[1] - FULL_Y[0]) / dy))
    y_full = yf0 + np.arange(ny_full) * dy
    i0 = int(round((y[0] - y_full[0]) / dy))          # computed window start in y
    out = []
    for s in slabs:
        base = np.full((ny_full, nx_full), plane_val, dtype=np.float32)
        base[i0:i0 + s.shape[0], j0:j0 + s.shape[1]] = s
        out.append(base)
    return out, x_full, y_full'''

assert src.count(old_helper) == 1, "pad_to_full body not matched verbatim"
src = src.replace(old_helper, new_helper)

# 3. update the call site: now returns 3 values and needs y passed in
old_call = 'slabs, xplot = pad_to_full(slabs, x, C0, mode, siz)'
new_call = 'slabs, xplot, yplot = pad_to_full(slabs, x, y, C0, mode, siz)'
assert src.count(old_call) == 1, "call site not unique"
src = src.replace(old_call, new_call)

# the animate call currently passes y; route the padded yplot through it
old_anim = 'animate(mode, sv, slabs, xplot, y, z_tr, fn, title_str=ttl, units_label=lab)'
new_anim = 'animate(mode, sv, slabs, xplot, yplot, z_tr, fn, title_str=ttl, units_label=lab)'
assert src.count(old_anim) == 1, "animate call not unique"
src = src.replace(old_anim, new_anim)
# and set yplot=y in the non-full-domain branch so the name always exists
old_xp = '        xplot = x\n'
new_xp = '        xplot = x; yplot = y\n'
assert src.count(old_xp) == 1, "xplot init not unique"
src = src.replace(old_xp, new_xp)

# 4. pin endpoints: x for all, y for xy plan views
old_xlim = '    ax.set_xlim(x[0], x[-1])'
new_xlim = '    ax.set_xlim(FULL_X[0], FULL_X[-1]) if FULL_DOMAIN else ax.set_xlim(x[0], x[-1])'
assert src.count(old_xlim) == 1, "set_xlim not unique"
src = src.replace(old_xlim, new_xlim)

old_ylim = '''    else:
        ax.set_ylim(V.min(), V.max())'''
new_ylim = '''    else:
        if FULL_DOMAIN:
            ax.set_ylim(FULL_Y[0], FULL_Y[-1])
        else:
            ax.set_ylim(V.min(), V.max())'''
assert src.count(old_ylim) == 1, "xy set_ylim not unique"
src = src.replace(old_ylim, new_ylim)

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_yfd")
open(SRC, "w").write(src)
print("OK, backup nutrient_tracer_fmt.py.bak_yfd")
