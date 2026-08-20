import ast, shutil, sys
SRC = "nutrient_tracer_fmt.py"
src = open(SRC).read()

def sub(old, new, what):
    global src
    n = src.count(old)
    if n != 1: sys.exit("ABORT: %s found %d times, expected 1" % (what, n))
    src = src.replace(old, new)

# 1. Nutrients -> Nitrate (display, colorbar, filenames)
sub('lab = "Nutrients (mmol m$^{-3}$)"; nut = "Nutrients"; fp = "nutrients"',
    'lab = "Nitrate (mmol m$^{-3}$)"; nut = "Nitrate"; fp = "nitrate"', 'labels')

# 2. X-Z y-label: z (m) only
sub('"z (m), 0 = surface\\n(AZMP 3-95 m compressed to 40 m)"',
    '"z (m)"', 'xz ylabel')

# 3. X-Z ticks including 0 (surface) and 40 (bed)
sub('        ax.set_ylim(V.max(), V.min())   # depth 0 (surface) top, 40 (seabed) bottom',
    '        ax.set_ylim(40.0, 0.0)          # 0 (surface) top, 40 (bed) bottom\n'
    '        ax.set_yticks(np.arange(0, 41, 5))', 'xz ylim/ticks')

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_nitrate")
open(SRC, "w").write(src)
print("OK, backup nutrient_tracer_fmt.py.bak_nitrate")
