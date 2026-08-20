W SPILL, WITH BUOYANCY
======================
See ../../ReadMe.txt for the velocity source, boundary conditions, cache
mechanics, and the meaning of the two limits.

PARAMETERS
  OIL_MASS_KG   10.0     0.0118 m3 at 850 kg/m3
  OIL_RADIUS    15.0 m
  OIL_OFFSET    24.0 m   spill centre to pile centre
  BEARING       W        spill centre (16, 60)
  VERTICAL_MIX  False    w zeroed at load
  TRACER_NZ     150      dz = 0.265 m

DERIVED
  Areal loading   14.15 g/m2
  Film thickness  16.6 um
  Surface cell    53.3 g/m3
  99.9% of the disc is inside the grid, 9.99 of 10.0 kg

Slick spans x = 1 to 31; the pile surface is at x = 36, so it starts
5 m clear and drifts in.

RESULT
  results/oil_W_buoyant.gif
  Colour range 0 to 15 g/m2, measured from the data.

Mass wanders about 20 percent below its initial value around mid-run and
returns close to 10 kg by the end. That is flow divergence in the binned
field, not a boundary problem.

The slick holds its peak concentration for roughly half the run, then the
flat-topped disc erodes from the edges inward. Footprint grows about 4x
while total mass is conserved, so the apparent fading is dilution.

BOUNDARY CODE
  Inlet fix present. No surface reflection - with w = 0 there is no
  vertical displacement, so there is nothing to reflect.
