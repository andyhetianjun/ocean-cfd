import ast, shutil, sys
SRC = "nutrient_tracer_fmt.py"
src = open(SRC).read()

def sub(old, new, what):
    global src
    n = src.count(old)
    if n != 1: sys.exit("ABORT: %s found %d times, expected 1" % (what, n))
    src = src.replace(old, new)

sub('def animate(mode, slice_val, slabs, x, y, z, filename, title_str="", fps=24, units_label=""):',
'''XZ_ASPECT     = 3.0
CBAR_LABELPAD = 3
CBAR_RANGES = {
    ("xz", None): (2.0,  8.0,  1.0),
    ("xy",  0):   (0.08, 0.24, 0.02),
    ("xy", 10):   (0.0,  1.2,  0.2),
    ("xy", 20):   (2.0,  6.0,  0.5),
    ("xy", 30):   (7.0, 10.0,  0.5),
    ("xy", 40):   (9.3, 10.0,  0.1),
}


def animate(mode, slice_val, slabs, x, y, z, filename, title_str="", fps=24, units_label=""):''',
    'config')

sub('''        figsize, xl, yl = (8, 4), "x (m)", "z (m)"
    else:
        H, V = np.meshgrid(x, y)
        figsize, xl, yl = ((14, 3.5) if FULL_DOMAIN else (9, 4)), "x (m)", "y (m)"
    lo = float(min(s.min() for s in slabs)); hi = float(max(s.max() for s in slabs))
    if hi - lo < 1e-9:
        lo, hi = lo - 0.05*abs(lo) - 1e-6, hi + 0.05*abs(hi) + 1e-6
    raw = (hi - lo) / 6.0
    mag = 10.0 ** np.floor(np.log10(raw))
    step = next(m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag)
    lo = np.floor(lo / step) * step
    hi = np.ceil(hi / step) * step
    cticks = np.round(np.arange(lo, hi + step * 0.5, step), 10)
    norm = Normalize(vmin=lo, vmax=hi)
    print(f"      colour range {lo:.3f} .. {hi:.3f}", flush=True)''',
'''        figsize, xl, yl = (14, 3.5), "x (m)", "z (m)"
        rkey = ("xz", None)
    else:
        H, V = np.meshgrid(x, y)
        figsize, xl, yl = ((14, 3.5) if FULL_DOMAIN else (9, 4)), "x (m)", "y (m)"
        rkey = ("xy", int(round((Z_SURFACE - slice_val) / 10.0)) * 10)
    rng = CBAR_RANGES.get(rkey)
    if rng is not None:
        lo, hi, step = rng
    else:
        lo = float(min(s.min() for s in slabs)); hi = float(max(s.max() for s in slabs))
        if hi - lo < 1e-9:
            lo, hi = lo - 0.05*abs(lo) - 1e-6, hi + 0.05*abs(hi) + 1e-6
        raw = (hi - lo) / 6.0
        mag = 10.0 ** np.floor(np.log10(raw))
        step = next(m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag)
        lo = np.floor(lo / step) * step
        hi = np.ceil(hi / step) * step
    lo = max(lo, 0.0)
    cticks = np.round(np.arange(lo, hi + step * 0.5, step), 10)
    norm = Normalize(vmin=lo, vmax=hi)
    _tot = sum(s.size for s in slabs)
    _clip = sum(int((s < lo).sum()) + int((s > hi).sum()) for s in slabs)
    _ext = "both" if _clip else "neither"
    print(f"      colour range {lo:.3f} .. {hi:.3f}   clipped {100.0*_clip/_tot:.2f}%", flush=True)''',
    'ranges')

sub('''    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"),
                 cax=cax, label=units_label, ticks=cticks)''',
'''    _cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"),
                       cax=cax, ticks=cticks, extend=_ext)
    _cb.set_label(units_label, labelpad=CBAR_LABELPAD)''',
    'colorbar')

sub('''        ax.set_yticks(np.arange(0, 41, 5))''',
'''        ax.set_yticks(np.arange(0, 41, 5))
        ax.set_aspect(XZ_ASPECT, adjustable="box")''',
    'aspect')

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_cbar")
open(SRC, "w").write(src)
print("OK, backup nutrient_tracer_fmt.py.bak_cbar")
