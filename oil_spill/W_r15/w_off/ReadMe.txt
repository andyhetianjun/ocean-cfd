OIL SPILL, SURFACE RELEASE, BEARING W
=====================================

Surface oil slick released upstream of the monopile and advected through
the LES wake. Passive tracer. No rise velocity, no weathering, no
evaporation, no dissolution.

CASE DEFINITION
---------------
Two limiting cases bracket the real behaviour.

  Case 1   VERTICAL_MIX = False
           w is zeroed. Oil stays in the surface layer.
           Units g/m2. Perfectly buoyant limit.

  Case 2   VERTICAL_MIX = True
           w is active, oil mixes downward.
           Units g/m3, because a surface slice reported in g/m2 stops
           meaning anything once material leaves that layer.
           Neutrally buoyant limit.

Real oil sits between the two. The model has no rise velocity, so neither
case is the answer on its own; they are the bounds.

Case 1 is complete. Case 2 has not been run.

PARAMETERS (current)
--------------------
  OIL_MASS_KG   10.0     0.0118 m3 at 850 kg/m3
  OIL_RADIUS    15.0 m
  OIL_OFFSET    24.0 m   spill centre to pile centre
  BEARING       W
  VERTICAL_MIX  False
  CYL_X, CYL_Y  40.0, 60.0
  CYL_R         4.0
  TRACER_NZ     150      dz = 0.265 m

Spill centre (16, 60). Slick spans x = 1 to 31. The pile surface is at
x = 36, so the slick starts 5 m clear and drifts in rather than beginning
in contact.

Derived: 14.15 g/m2 areal loading, 53.3 g/m3 in the surface cell,
16.6 um film thickness. The film thickness sits between the two cases in
Yongsheng's original table and does not enter the model, which works in
mass per unit area.

PARAMETER HISTORY
-----------------
Yongsheng's original table gave 0.1 m3, 85 kg, and either r = 56.4 m at
10 um or r = 17.8 m at 100 um. A 56.4 m radius does not fit a 120 m wide
tank. Revised to r = 20 m at 24 m offset, then to r = 15 m and 10 kg after
the pile position was corrected.

The pile correction matters more here than for nitrate. Scripts inherited
CYL_X = 46.0 from older handoff notes; the mesh says x = 40. Because the
spill is placed relative to the pile, the wrong value put the first run
18 m from the pile instead of 24, with part of the slick initialised inside
the pile where there is no flow. That was the "oil stuck to the structure"
behaviour.

VELOCITY
--------
velocity_full_025, 0.25 m binning, grid 2400 x 480 x 84.
Timesteps 1049-1400 s, 352 steps.
Same source and time block as the nitrate work. The oil spill was
originally built against the superseded 399-595 s block; it was moved to
the full-domain extraction because the W spill sits upstream at x = 16 and
needs coverage from x = 1.

Resolution limit is inherited: the mesh is refined to x approx 290 m and
coarsens beyond, so downstream of x approx 350 m the figures show
large-scale displacement rather than resolved filaments.

OUTPUT
------
  results/oil_W_xy_z0_surface_aspect1p5.gif     current, for Yongsheng
  results/oil_W_xy_z0_surface_aspect_equal.gif  same run, true 1:1

Colour range 0.000 to 15.000, measured off the data rather than fixed. The
script computes lo/hi from the slabs and snaps to a round step, so the range
tracks the initial loading automatically. Earlier runs at 85 kg produced a
0-80 range; that is not a hardcoded value and needed no editing.

ASPECT
------
The plan view is now ax.set_aspect(1.5, adjustable="box"), matching the
nitrate plan views at Yongsheng's request.

Worth stating whenever this figure goes out: at 1.5 the slick and the pile
both render as ellipses, and cross-stream distance is stretched 1.5x
relative to downstream. The slick is a hard-edged circle of known radius,
so this costs more here than it does for nitrate, where there is no sharp
geometry to distort. The equal-aspect version is kept alongside for that
reason. It is a display choice and does not change the spill geometry.

RUNNING
-------
Run from this folder, not from scripts/.

  /home/andyhe/tracer_env/bin/python -u scripts/oil_spill.py > render_oil.log 2>&1 &

Check free -g first; this needs 20-30 GB on a shared machine.

Runtime: the full advection took 25863 s, about 7.2 hours. Note this is
well above the 3 hours quoted in earlier notes.

Cached re-encodes take about 7 minutes. The cache is
data/slabs_cache/arr_0.npy, a single stacked array of about 1.6 GB, plus
sig.txt. It is not per-frame npz files.

Signature format in sig.txt:
  velocity dir basename | TRACER_NZ | first step | last step | n steps |
  IC hash | VERTICAL_MIX | slice spec

Aspect and colour range are not in the signature, so both are display-only
and re-encode off cache.

OUTPUT PATH
-----------
The animation save previously took a bare filename, so output landed in
this folder rather than in results/ and had to be moved by hand. Fixed:
the save now uses os.path.join("results", filename).

Consequence: consecutive runs overwrite in place. BEARING is in the
filename so different bearings are safe, but case 2 at the same bearing
would clobber case 1. Move or rename before running case 2.

BEARINGS
--------
Eight bearings were planned at 45 degree increments. Each is a separate
initial condition, so the slab cache cannot be shared between them.

At 7.2 hours each that is about 58 hours serial, not the 24 hours in
earlier notes. Worth confirming with Yongsheng how many he actually wants
before running them all. E, NE and SE start downstream of the pile and
interact with the wake much less; N and S may clip only the wake edge.

TODO
----
- Confirm how many bearings are wanted.
- Case 2 (VERTICAL_MIX = True) not yet run.
