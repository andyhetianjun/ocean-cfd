import ast, shutil, sys
SRC = "nutrient_tracer_fmt.py"
src = open(SRC).read()

def sub(old, new, what):
    global src
    n = src.count(old)
    if n != 1: sys.exit("ABORT: %s found %d times, expected 1" % (what, n))
    src = src.replace(old, new)

sub('at 0 m depth (surface)',     'at z = 0 m (surface)',     'd00 title')
sub('at 10 m depth (subsurface)', 'at z = 10 m (subsurface)', 'd10 title')
sub('at 20 m depth (middle)',     'at z = 20 m (middle)',     'd20 title')
sub('at 30 m depth (sub-bottom)', 'at z = 30 m (sub-bottom)', 'd30 title')
sub('at 40 m depth (bottom)',     'at z = 40 m (bottom)',     'd40 title')
sub('"depth below surface (m)',   '"z (m), 0 = surface',      'xz ylabel')

sub('    norm = Normalize(vmin=lo, vmax=hi)',
'''    raw = (hi - lo) / 6.0
    mag = 10.0 ** np.floor(np.log10(raw))
    step = next(m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag)
    lo = np.floor(lo / step) * step
    hi = np.ceil(hi / step) * step
    cticks = np.round(np.arange(lo, hi + step * 0.5, step), 10)
    norm = Normalize(vmin=lo, vmax=hi)''', 'norm block')

sub('                 cax=cax, label=units_label)',
    '                 cax=cax, label=units_label, ticks=cticks)', 'colorbar ticks')

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_ticks")
open(SRC, "w").write(src)
print("OK, backup nutrient_tracer_fmt.py.bak_ticks")
