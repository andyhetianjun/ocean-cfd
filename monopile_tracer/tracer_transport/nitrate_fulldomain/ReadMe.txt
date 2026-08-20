NITRATE TRACER, FULL DOMAIN
===========================

Advection of an AZMP nitrate profile through the LES wake of an 8 m monopile.
Passive tracer, advection only, no diffusion and no biology.

SOURCE
------
LES velocity from
  /shared_folder/yongxing/OpenFOAM/simulationCases/flow_past_cylinder/
      realistic_cases/domain_D8X600Y120/uniform
Domain 600 x 120 x 40 m. Monopile D = 8 m centred at x = 40, y = 60,
spanning the full depth. Freestream approx 0.5 m/s, vertically uniform,
slip bed.

Extraction: velocity_full_025, 0.25 m horizontal binning over the full
domain, grid 2400 x 480 x 84. Timesteps 1049-1400 s, 352 steps.
Empty bins are filled by nearest neighbour; approx 60 percent of cells at
0.25 m are filled this way.

RESOLUTION LIMIT
----------------
The LES mesh is refined to x approx 290 m and coarsens to 1-3 m beyond.
Past x approx 350 m the figures show large-scale displacement of the tracer,
not resolved turbulent filaments. Confirmed as expected by Yongxing.

INITIAL CONDITION
-----------------
AZMP discrete-bottle nitrate, nine binned depths, stretched over the 40 m
domain. Horizontally uniform at t = 0.
Interpolation is np.interp followed by np.maximum.accumulate, which
guarantees a monotonic profile.
An earlier version used UnivariateSpline(s=0.0), which overshot between
knots and put the true minimum 4.8 m below the surface and the true maximum
2.4 m above the bed. That made the z = 0 and z = 40 planes show values
outside their own backgrounds. Advection was never at fault; the
semi-Lagrangian trilinear scheme is bounded by construction.

Tracer grid: 100 vertical levels.

NUTRICLINE
----------
Peak gradient        17.4 m
Base                 33.8 m   (20 percent of peak gradient; a chosen
                               threshold, not a standard convention)
Top                  12 m     (Yongsheng's read by eye)
Three plan views are rendered at 12, 17 and 20 m because there are three
defensible answers to where the nutricline sits.

OUTPUT
------
Thirteen animations in results/, 24 fps, absolute simulation timestamps.

  X-Z sections (5), at y =
     56, 58, 60, 62, 64          colour range 2-8 for all five

  X-Y plan views (8), at z =
     0    colour range 0.08-0.24
     10   colour range 0-1.2
     12   colour range 0-2.0
     17   colour range 1.5-5.0
     20   colour range 2-6
     30   colour range 7-10
     34   colour range 7-10  (inherited, see note below)
     40   colour range 8-10

The section range of 2-8 clips roughly 65 percent of the field. This is
deliberate: the flat caps at the top and bottom of the profile are cut so
that the mixing layer fills the colour bar.

Also in results/:
  initial_profile.png        the AZMP profile as applied
  full_domain_velocity.png   velocity field over the full domain

DISPLAY SETTINGS
----------------
XZ_ASPECT = 4.5   vertical exaggeration on sections
XY_ASPECT = 1.5   on plan views

Both constants sit near the top of scripts/nutrient_tracer_fmt.py, just
above CBAR_RANGES, and are applied in animate() via
ax.set_aspect(N, adjustable="box").

History of those values: the sections originally had no set_aspect call at
all, which measured out at approx 8.5x. Yongsheng asked for a third of that,
giving 3.0; he then found 3.0 too flat and asked to move 25 percent back
toward the original, giving 3 + 0.25 * (8.54 - 3) = 4.4, rounded to 4.5.
The plan views were originally set_aspect("equal"), true 1:1, which on a
600 x 120 domain is a 5:1 strip that reads badly on a slide. 1.5 was chosen
to match the sections. Tradeoff: the monopile draws as an ellipse rather
than a circle, and cross-stream distance is stretched 1.5x relative to
downstream.

Pile drawn in black, no white ring.
Timestamp bottom-left on plan views, bottom-right on sections.

RUNNING
-------
Run from this folder, not from scripts/, because paths are relative.

  /home/andyhe/tracer_env/bin/python -u scripts/nutrient_tracer_fmt.py

Check free -g first. The 2400 x 480 grid needs 20-30 GB and the machine is
shared; three runs have been lost to the OOM killer because of other users'
jobs. The killer picks the largest process, which is usually this one.

SLAB CACHE
----------
evolve() caches advected slices to data/slabs_cache/. Display-only changes
(colour ranges, aspect, labels, fps) then re-encode instead of re-advecting.
The cache signature covers: velocity directory basename, TRACER_NZ, timestep
range, slice list, initial-condition hash, and the vertical-mixing flag.
Changing any of those forces a full re-advection. Aspect and colour range
are not in the signature, which is why they are cheap to change.

COLOUR RANGE LOOKUP
-------------------
Ranges live in CBAR_RANGES in scripts/nutrient_tracer_fmt.py, keyed by
(mode, depth). There is no explicit entry for z = 34. The lookup at lines
177-178 falls back to the nearest existing key, so 34 takes z = 30's range
of 7-10. That is the intended range, but it is reached by fallback rather
than by an explicit entry.

The fallback is silent. A new plane added at a depth with no key will
inherit its nearest neighbour's range instead of auto-measuring, and
nothing in the output says so. Check CBAR_RANGES before adding planes.

FILES IN results/
-----------------
  nitrate_xy_d00_surface.gif
  nitrate_xy_d10_subsurface.gif
  nitrate_xy_d12_nutriclinetop.gif
  nitrate_xy_d17_nutriclinepeak.gif
  nitrate_xy_d20_middle.gif
  nitrate_xy_d30_subbottom.gif
  nitrate_xy_d34_nutriclinebase.gif
  nitrate_xy_d40_bottom.gif
  nitrate_xz_y56.gif
  nitrate_xz_y58.gif
  nitrate_xz_y60.gif
  nitrate_xz_y62.gif
  nitrate_xz_y64.gif
  initial_profile.png
  full_domain_velocity.png

HISTORY
-------
Pile position was previously CYL_X = 46.0, inherited from older notes; the
mesh says 40. Nitrate concentrations were unaffected, only the drawn marker
moved. Any figures in archive/ predating this show the marker at 46.
