#!/usr/bin/env python3
"""
Plot the tracer initial condition: the vertical nitrate profile as mapped
onto the 40 m model domain, with the nutricline marked.

Run from the figure directory:
    python scripts/plot_initial_profile.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# pull the profile-building code out of the tracer script so this figure
# can never drift from what the model actually uses
_src = open("scripts/nutrient_tracer_fmt.py").read().split("def load_grid")[0]
exec(_src)

OUT = "initial_profile.png"

pdp, pvl = azmp_profile(AZMP_CSV, NUTRIENT)
import glob; _zv = np.load(sorted(glob.glob(VEL_DIR + "/vel_*.npz"))[0])["z"]; z = np.linspace(_zv[0], _zv[-1], TRACER_NZ)
C0, depth_equiv = build_initial_field(z, pdp, pvl, 1, 1)

c = C0[:, 0, 0][::-1]                 # surface -> bed
d = (Z_SURFACE - z)[::-1]
g = np.gradient(c, d)

pk = int(np.argmax(g))
thr = 0.2 * g[pk]
below = np.where(g[pk:] < thr)[0]
base = float(d[pk + below[0]]) if below.size else float(d[-1])

print("surface %.3f   bed %.3f" % (c[0], c[-1]))
print("peak gradient %.3f at %.1f m" % (g[pk], d[pk]))
print("nutricline base %.1f m" % base)

fig, ax = plt.subplots(figsize=(5.5, 6), constrained_layout=True)

ax.axhspan(0, base, color="#2a78d6", alpha=0.06, zorder=0)
ax.plot(c, d, color="#2a78d6", lw=2.2, zorder=3)
ax.axhline(12.0, color="#eb6834", lw=1.2, ls=":", zorder=2)
ax.text(9.9, 12.0 - 0.6, "top of nutricline, 12.0 m",
        ha="right", va="bottom", fontsize=9, color="#eb6834")
ax.axhline(d[pk], color="#eb6834", lw=1.2, ls="--", zorder=2)
ax.text(9.9, d[pk] - 0.6, "steepest gradient, %.1f m" % d[pk],
        ha="right", va="bottom", fontsize=9, color="#eb6834")
ax.axhline(base, color="#eb6834", lw=1.2, ls=":", zorder=2)
ax.text(9.9, base - 0.6, "base of nutricline, %.1f m" % base,
        ha="right", va="bottom", fontsize=9, color="#eb6834")
ax.set_xlim(0, 10.4)
ax.set_ylim(40, 0)
ax.set_xlabel("nitrate (mmol m$^{-3}$)")
ax.set_ylabel("z (m), 0 = surface")
ax.set_yticks(np.arange(0, 41, 5))
ax.grid(alpha=0.25, lw=0.6)

axr = ax.twinx()
axr.set_ylim(40, 0)
ticks = np.arange(0, 41, 5)
axr.set_yticks(ticks)
axr.set_yticklabels(["%.0f" % (AZMP_TOP + (t / 40.0) * (AZMP_BOT - AZMP_TOP))
                     for t in ticks], fontsize=9)
axr.set_ylabel("corresponding AZMP depth (m)", fontsize=10)

ax.set_title("Nitrate initial condition\nAZMP station HL5 (Jun-Sep), stretched "
             "from %g-%g m into the %g m domain"
             % (AZMP_TOP, AZMP_BOT, Z_SURFACE - Z_BOTTOM),
             fontsize=10.5, fontweight="bold")

fig.savefig(OUT, dpi=150)
print("saved", OUT)
