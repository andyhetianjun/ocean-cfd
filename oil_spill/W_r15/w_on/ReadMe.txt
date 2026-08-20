W SPILL, WITH NO BUOYANCY
=========================
See ../../ReadMe.txt for the velocity source, boundary conditions, cache
mechanics, and the meaning of the two limits.

Identical to ../w_off except for VERTICAL_MIX, which is the point - the
pair is meant to be compared directly.

PARAMETERS
  OIL_MASS_KG   10.0     0.0118 m3 at 850 kg/m3
  OIL_RADIUS    15.0 m
  OIL_OFFSET    24.0 m
  BEARING       W        spill centre (16, 60)
  VERTICAL_MIX  True     w active
  TRACER_NZ     150      dz = 0.265 m

DERIVED
  Areal loading   14.15 g/m2 at t = 0
  Surface cell    53.3 g/m3
  99.9% of the disc is inside the grid, 9.99 of 10.0 kg

RESULT
  results/oil_W_nonbuoyant.gif
  Colour range 0 to 60 g/m3.

Units are g/m3 rather than g/m2 because oil leaves the surface layer, so
an areal density on a single slice is meaningless.

The surface slice fades through the run as oil mixes downward. Nothing
returns it - there is no rise velocity in the model - so this is the
maximum possible surface clearing.

Mass runs 10.04 to 10.64 kg, about +6 percent.


BOUNDARY CODE
  Inlet fix present. Surface reflection present.
