import ast, shutil, sys
SRC = "oil_spill.py"
src = open(SRC).read()

LAB_OLD = '    lab = "Nitrate (mmol m$^{-3}$)"; nut = "Nitrate"; fp = "nitrate"'
RET_OLD = '''    return (np.nan_to_num(g[0]).astype(np.float32),
            np.nan_to_num(g[1]).astype(np.float32),
            np.nan_to_num(g[2]).astype(np.float32))'''
IC_OLD = '''    pdp, pvl = azmp_profile(AZMP_CSV, NUTRIENT)
    C0, depth_equiv = build_initial_field(z_tr, pdp, pvl, ny, nx)
    print(f"IC {NUTRIENT}: {C0.min():.3f} .. {C0.max():.3f}", flush=True)'''
SIG_OLD = '''    sig = "|".join([VEL_DIR, str(TRACER_NZ), TIMESTEPS[0], TIMESTEPS[-1],
                    str(len(TIMESTEPS))] + [f"{m}:{sv}" for m, sv, _, _ in jobs])'''

anchors = {
    "VEL_DIR":   'VEL_DIR    = "velocity_resolved"',
    "TRACER_NZ": 'TRACER_NZ = 150',
    "load_grid": 'def load_grid():',
    "uvw_return": RET_OLD,
    "nitrate_IC": IC_OLD,
    "lab_line":  LAB_OLD,
    "cache_sig": SIG_OLD,
}
missing = {k: src.count(v) for k, v in anchors.items() if src.count(v) != 1}
if missing:
    print("ABORT, anchors not found exactly once:")
    for k, n in missing.items():
        print("  %-12s count=%d" % (k, n))
    sys.exit(1)

def sub(old, new):
    global src
    src = src.replace(old, new)

sub('VEL_DIR    = "velocity_resolved"',
    'VEL_DIR    = "/shared_folder/andyhe/project/tracer_transport/velocity_resolved"')

sub('TRACER_NZ = 150',
'''TRACER_NZ = 150

# --- oil spill config -------------------------------------------------
OIL_MASS_KG   = 85.0        # 0.1 m3 at 850 kg/m3
OIL_RADIUS    = 20.0        # m
OIL_OFFSET    = 24.0        # m, spill centre to pile centre
BEARING       = "W"         # N NE E SE S SW W NW
VERTICAL_MIX  = False       # case 1: w = 0.  case 2: True
_BEAR = {"N": (0, 1), "NE": (0.7071, 0.7071), "E": (1, 0), "SE": (0.7071, -0.7071),
         "S": (0, -1), "SW": (-0.7071, -0.7071), "W": (-1, 0), "NW": (-0.7071, 0.7071)}
OIL_CX = CYL_X + _BEAR[BEARING][0] * OIL_OFFSET
OIL_CY = CYL_Y + _BEAR[BEARING][1] * OIL_OFFSET
# ----------------------------------------------------------------------''')

sub('def load_grid():',
'''def build_oil_field(x, y, z_tr):
    """Disc of oil in the surface cell. Returns (C in g/m3, diagnostics)."""
    nz, ny, nx = len(z_tr), len(y), len(x)
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0]); dz = float(z_tr[1] - z_tr[0])
    X, Y = np.meshgrid(x, y)
    inside = ((X - OIL_CX)**2 + (Y - OIL_CY)**2) <= OIL_RADIUS**2
    disc_area = np.pi * OIL_RADIUS**2
    conc = (OIL_MASS_KG * 1000.0) / (disc_area * dz)      # g/m3
    C = np.zeros((nz, ny, nx), dtype=np.float32)
    C[-1][inside] = conc                                   # z_tr[-1] is the surface
    grid_area = float(inside.sum()) * dx * dy
    return C, {"conc_gm3": conc, "areal_gm2": conc * dz,
               "frac_in_grid": grid_area / disc_area,
               "mass_in_grid_kg": conc * grid_area * dz / 1000.0}


def load_grid():''')

sub(RET_OLD,
'''    w = np.nan_to_num(g[2]).astype(np.float32)
    if not VERTICAL_MIX:
        w = np.zeros_like(w)
    return (np.nan_to_num(g[0]).astype(np.float32),
            np.nan_to_num(g[1]).astype(np.float32),
            w)''')

sub(IC_OLD,
'''    C0, oil = build_oil_field(x, y, z_tr)
    print(f"oil: r={OIL_RADIUS} m at ({OIL_CX:.1f}, {OIL_CY:.1f}) bearing {BEARING}", flush=True)
    print(f"     {oil['conc_gm3']:.1f} g/m3 in surface cell = {oil['areal_gm2']:.2f} g/m2", flush=True)
    print(f"     {oil['frac_in_grid']*100:.1f}% of disc inside grid, "
          f"{oil['mass_in_grid_kg']:.2f} of {OIL_MASS_KG} kg", flush=True)
    print(f"     vertical mixing: {'ON (case 2)' if VERTICAL_MIX else 'OFF (case 1, w=0)'}", flush=True)''')

sub(LAB_OLD, '    lab = "Oil (g m$^{-2}$)"; nut = "Oil"; fp = f"oil_{BEARING}"')

sub(SIG_OLD,
'''    import hashlib
    ich = hashlib.md5(C0.tobytes()).hexdigest()[:12]
    sig = "|".join([VEL_DIR, str(TRACER_NZ), TIMESTEPS[0], TIMESTEPS[-1],
                    str(len(TIMESTEPS)), ich, str(VERTICAL_MIX)]
                   + [f"{m}:{sv}" for m, sv, _, _ in jobs])''')

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_oil")
open(SRC, "w").write(src)
print("OK")
