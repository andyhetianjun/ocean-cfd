OIL SPILL AT THE MONOPILE
=========================

Four cases. Surface oil advected through the LES wake of an 8 m monopile.
Passive tracer: no rise velocity, no weathering, no evaporation, no
dissolution, no diffusion term.

  W_r15/w_off        disc spill west of the pile, with buoyancy
  W_r15/w_on         same disc, with no buoyancy
  fulldomain/w_off   oil covering the whole surface, with buoyancy
  fulldomain/w_on    same, with no buoyancy

Each folder has its own ReadMe with parameters and results.

THE TWO LIMITS
--------------
The model has no buoyancy force and no rise velocity. The two cases impose
the two extremes rather than predicting anything in between.

  w_off  VERTICAL_MIX = False. Vertical velocity is zeroed at load, so oil
         stays in the surface layer. Reported in g/m2. This is the
         perfectly buoyant limit.

  w_on   VERTICAL_MIX = True. Vertical velocity is active and oil mixes
         downward freely. Reported in g/m3, because a surface slice in
         g/m2 stops meaning anything once material leaves that layer.
         This is the neutrally buoyant limit.

Real oil sits between the two. The model cannot say where.

Note the phrasing: oil stays at the surface in w_off because w was set to
zero, not because the model knows oil is buoyant. It is an imposed
assumption, not a result.

VELOCITY SOURCE
---------------
  /shared_folder/yongxing/OpenFOAM/simulationCases/flow_past_cylinder/
      realistic_cases/domain_D8X600Y120/uniform

Extracted to a regular grid by scripts/extract_full_025.py, giving
data/velocity_full_025: 2400 x 480 x 84 at 0.25 m horizontal spacing,
timesteps 1049-1400 s (352 steps).

Domain 600 x 120 x 40 m. Monopile D = 8 m centred at x = 40, y = 60,
spanning the full depth. Freestream approx 0.5 m/s, vertically uniform,
slip bed. Side walls at y = 0 and y = 120 are slip or symmetry: measured
cross-stream velocity there is about 2e-5 m/s against 0.067 m/s
mid-domain, so nothing crosses them.

Grid z runs -19.76 to +19.76; the surface is the last index.

RESOLUTION LIMIT
----------------
The LES mesh is refined to x approx 290 m and coarsens to 1-3 m beyond.
Past roughly x = 350 m the figures show large-scale displacement rather
than resolved turbulent filaments. Confirmed as expected by Yongxing.
About 60 percent of cells in the extraction are filled by nearest
neighbour.

BOUNDARY CONDITIONS
-------------------
The semi-Lagrangian step traces each cell back to where its water came
from one timestep earlier. Departure points that fall outside the domain
need explicit handling, and each face is treated differently.

Inlet at x = 0.
  Disc cases: tracer entering through the inlet is set to zero, since a
  finite spill has nothing upstream of it.
  Full-domain cases: the inlet holds its initial value, so the sheet is
  continuous. Without this the sheet drains from x = 0 at 0.5 m/s and the
  pile, at x = 40, sits inside the drained region within about 80
  seconds.

Side walls at y = 0 and y = 120.
  Clipped. These are slip or symmetry boundaries in the LES, so nothing
  crosses them and clipping approximates a reflecting condition well
  enough at this resolution.

Free surface and bed.
  Reflected, not clipped. Where vertical velocity is downward, the
  departure point lies above the free surface; the parcel actually came
  from just below it. Clipping would return the point to the surface cell
  itself, which duplicates tracer wherever there is downwelling - about
  57 percent of surface cells, with vertical displacements up to 0.10 m
  against a cell height of 0.265 m. Reflection is the correct treatment
  for a no-flux boundary and conserves mass.

  This matters only when VERTICAL_MIX is on. The w_off cases have no
  vertical displacement, so their scripts have no reflection.

MASS CONSERVATION
-----------------
Semi-Lagrangian trilinear interpolation is bounded but not exactly
conservative, and the binned velocity field is not perfectly
divergence-free, so some drift is expected even with correct boundaries.

  W with buoyancy         wanders about 20 percent and returns
  W with no buoyancy      +6 percent over the run
  Full domain, buoyancy   flat to 0.1 percent
  Full domain, no buoy.   +4 percent over the run


The nitrate animations in ../tracer_transport/ use the same routine with
vertical velocity active and drift 1.6 percent. The effect is much
smaller there because the nitrate profile is spread through the water
column rather than concentrated in the surface cell where the clipping
acted. Those figures have not been re-run.

SLAB CACHE
----------
evolve() caches advected slices to data/slabs_cache/ as a single stacked
arr_0.npy plus sig.txt. Display-only changes - colour range, aspect,
titles, frame rate - then re-encode in about 8 minutes instead of
re-advecting for 7 hours.

Signature format:
  velocity dir | TRACER_NZ | first step | last step | n steps |
  IC hash | VERTICAL_MIX | slice spec

Anything not in that list is free to change. Anything in it forces a full
re-advection.

RUNNING
-------
Run from the case folder, not from scripts/, because paths are relative.

  cd W_r15/w_off
  setsid nohup /home/andyhe/tracer_env/bin/python -u scripts/oil_spill.py \
      > render.log 2>&1 < /dev/null &

Check free -g first. Each run needs 20-30 GB at the startup peak and the
machine is shared; several runs have been lost to the OOM killer, which
picks the largest process. Full advection is about 7 hours.

Note that setsid makes the shell print "Done" immediately. That means
detached, not finished. Use pgrep -af oil_spill.py.

DISPLAY
-------
Frame rate 12 fps (GIF stores the delay in hundredths of a second, so it
is written as 80 ms, giving 12.5). Plan view aspect 1.5, matching the
nitrate figures - the slick and pile therefore draw as ellipses, with
cross-stream distance stretched 1.5x relative to downstream. Pile drawn
in black. Colour ranges are measured from the data rather than fixed.

Output goes to results/ automatically.

KNOWN LIMITATIONS
-----------------
- Four near-identical copies of oil_spill.py, one per case, differing in
  a handful of constants. Every fix has to be applied four times, which
  is how the boundary problems above took as long as they did to track
  down. Worth collapsing into one script plus four config files.
- Concentration values are approximate at the few-percent level.
- Nothing constrains where between the two limits real oil sits.
