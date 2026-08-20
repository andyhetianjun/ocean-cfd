MONOPILE TRACER
===============

Passive tracer transport through the wake of an offshore monopile, using
archived LES velocity. The LES is not part of this work; everything here
is post-processing of it.

Four studies share the same velocity field, and the first three share the
same advection code:

  tracer_transport/nitrate_fulldomain/    a nutrient profile through the wake
  tracer_transport/vertical_diffusivity/  flux-gradient estimate of Kz
  oil_spill/                              surface oil, four cases
  vorticity/                              vorticity and speed animations

Each has its own ReadMe with parameters, results, and how to run it.


THE VELOCITY FIELD
------------------
Large-eddy simulation by Yongxing Ma:

  /shared_folder/yongxing/OpenFOAM/simulationCases/flow_past_cylinder/
      realistic_cases/domain_D8X600Y120/uniform

  Domain        600 x 120 x 40 m
  Monopile      D = 8 m at (40, 60), spanning the full depth
  Freestream    approx 0.5 m/s, vertically uniform
  Bed           slip, so no bottom boundary layer
  Sides         slip or symmetry at y = 0 and y = 120. Measured
                cross-stream velocity there is about 2e-5 m/s against
                0.067 m/s mid-domain, so nothing crosses them.
  Storage       decomposed across 256 processors, about 1.87 TB

The raw output is on an unstructured mesh. Everything here works from a
regular grid binned out of it by the scripts in scripts/.


EXTRACTION
----------
All five scripts do the same job with different windows and resolutions.
They read the decomposed OpenFOAM output in parallel, bin cell values onto
a regular grid, and write one npz per timestep. All reuse
velocity_450_cache.npz, which holds the cell centres - the mesh is static,
so they are read once and reused for every timestep.

  extract_full_025.py   0 to 600 m, 2400 x 480 x 84 at 0.25 m.
                        This produced velocity_full_025, which everything
                        current uses.
  extract_full.py       same extent at 0.5 m, 1200 x 240 x 84. Superseded.
  extract_new.py        0 to 290 m, y 15 to 115. An earlier window that
                        stopped at the refined-mesh limit.
  extract_resolved.py   20 to 290 m. Excluded the inflow region.
  extract_resolved_x0.py  as above but starting at x = 0, added when the
                        oil spill needed coverage upstream of the pile.

The progression runs from narrow and coarse to full-domain and fine. The
later scripts exist because each study needed more of the domain than the
one before it.

Empty bins are filled by nearest neighbour. At 0.25 m about 60 percent of
cells are filled this way, which is worth remembering for any quantity
that involves a spatial derivative - vorticity in particular is noisier
than one computed on the native mesh would be.


RESOLUTION LIMIT
----------------
The LES mesh is refined out to x of about 290 m and coarsens to 1-3 m
beyond. Past roughly x = 350 m the figures show large-scale displacement
of tracer rather than resolved turbulent filaments. Confirmed as expected
by Yongxing.

In pile diameters that is a usable range of about x/D = 35, which is a
reasonable near-to-intermediate wake but not a far-wake study.


THE ADVECTION SCHEME
--------------------
Eulerian, continuous, semi-Lagrangian with trilinear interpolation. No
diffusion term - background diffusivity is set to zero, so any mixing
that appears is resolved transport rather than something prescribed.

The scheme is bounded by construction, so a tracer cannot overshoot its
own initial range. It is not exactly conservative, and the binned
velocity field is not perfectly divergence-free, so total mass drifts by
a few percent over a 352-step run.

Boundary treatment differs by face and matters more than it looks:

  Inlet     zero inflow for a finite release, or holds its initial value
            for a continuous one. Which is correct depends on the case.
  Sides     clipped. They are no-flux in the LES, so this is adequate.
  Surface   reflected, not clipped. A parcel whose departure point lies
  and bed   above the free surface actually came from just below it.
            Clipping returns it to the surface cell, which duplicates
            tracer wherever there is downwelling. Reflection is the
            correct treatment for a no-flux boundary and conserves mass.

The surface treatment only matters when vertical velocity is active. It
matters a lot when it does: about 57 percent of surface cells have
downward w, with displacements up to 0.10 m against a cell height of
0.265 m.


WHAT THE STUDIES FOUND
----------------------

Oil spill.
  Two limits bracket the behaviour of a surface slick, because the model
  has no rise velocity and cannot say where between them real oil sits.

  With buoyancy imposed - vertical velocity zeroed - a slick released
  upstream stays at the surface and is stretched and folded by the wake.
  Its footprint grows substantially over the run while total mass is
  conserved, so the apparent fading is dilution rather than loss.

  With no buoyancy, the pile clears essentially all the oil out of the
  surface layer in the near wake. Surface concentrations in the wake core
  fall close to zero. Nothing returns it, so this is the maximum possible
  surface clearing.

  A slick covering the entire surface shows no signature at all under
  horizontal advection alone. That is not a null result to be explained
  away: translating a horizontally uniform field leaves it uniform, so
  there is no edge to distort and no gradient to stretch. The same case
  with vertical mixing active produces the clearest figure in the set -
  the wake reads as a depletion track through an otherwise untouched
  sheet, widening downstream and breaking into discrete vortex
  structures.

  Taken together: a surface slick passing a monopile is affected mainly
  through vertical mixing, not horizontal stirring. How much depends
  entirely on the oil's buoyancy, which this model does not represent.

Vertical diffusivity.
  Method implemented following sections 5.1 to 5.4 of Yongsheng's
  document: Reynolds decomposition of vertical velocity and tracer,
  layer-averaged turbulent flux, then a least-squares flux-gradient fit.

  A short test run gives Kz of order 1e-3 m2/s, uniform to within a
  factor of two over the water column. That figure is provisional - it
  comes from 30 of 352 timesteps and the full runs are still going.

  Two things are already clear from the test. Levels near the surface and
  bed carry no information, because the cosine initial profile has zero
  gradient there, so they have to be masked rather than reported. And the
  choice of averaging operator matters: a time mean and a layer mean give
  diffusivities differing by a factor of about 1.7. Which one section 5.1
  intends is an open question with Yongsheng.

Nitrate.
  Thirteen animations of an AZMP nitrate profile advected through the
  wake, produced to Yongsheng's specification - he chose the nutricline
  depths and the plane list. Interpretation of what they show is his;
  this repository holds the method and the output.

Vorticity and speed.
  Two-panel animations built for Gloria Wang at Yongxing's request, using
  the same extraction. Plan view at mid-depth over a streamwise section
  along the pile centreline.


PRACTICAL NOTES
---------------
Python environment: /home/andyhe/tracer_env/bin/python

Scripts are run from their case folder, not from scripts/, because the
paths inside them are relative.

Advection runs need 20-30 GB at their startup peak and take about seven
hours for 352 steps on a 2400 x 480 grid. The machine is shared; check
free -g before launching, since the OOM killer picks the largest process
and several runs have been lost that way.

Run detached so a job survives the terminal closing:

  setsid nohup /home/andyhe/tracer_env/bin/python -u scripts/foo.py \
      > render.log 2>&1 < /dev/null &

The shell prints "Done" immediately when you do this. That means
detached, not finished. Use pgrep to see what is actually running.

Print a running diagnostic every N steps - a total, a maximum, anything.
Three separate boundary-condition faults in this code were caught purely
because a printed mass total was climbing when it should have been flat.
None of them would have been visible in the figures.

Velocity data, advection caches and rendered animations are excluded from
version control for size. The scripts regenerate them.
