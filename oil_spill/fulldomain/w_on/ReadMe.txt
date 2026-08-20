FULL-DOMAIN SLICK, WITH NO BUOYANCY
===================================
See ../../ReadMe.txt for the velocity source, boundary conditions, cache
mechanics, and the meaning of the two limits.

This is the most informative of the four figures.

PARAMETERS
  OIL_MASS_KG      573.0    over 72000 m2
  OIL_FULL_DOMAIN  True
  VERTICAL_MIX     True     w active
  TRACER_NZ        150      dz = 0.265 m

  OIL_RADIUS, OIL_OFFSET and BEARING are unused. BEARING still appears in
  the log header - ignore it there.

DERIVED
  Areal loading   7.96 g/m2 at t = 0
  Surface cell    30.0 g/m3

RESULT
  results/oil_fulldomain_nonbuoyant.gif
  Colour range 0 to 40 g/m3.

The wake appears as depletion rather than concentration. The pile mixes
oil downward out of the surface layer, leaving a dark track through the
uniform background that widens downstream and breaks into discrete vortex
structures. Concentrations in the wake core fall close to zero, so in
this limit the pile clears essentially all the oil from the surface in
the near wake.

Reading the figure: the uniform background is untouched oil at its
initial value. Dark is oil that has been removed from the surface, not
oil that has not arrived. This is the opposite convention to the disc
figures.

Structure past about x = 350 m is coarser and blobbier. That is the mesh
coarsening from 0.25 m to 1-3 m, not a change in the physics.

Mass runs 571 to 594 kg, about +4 percent. Some increase is expected here
by design, since the inlet supplies fresh oil continuously while only
what reaches x = 600 leaves.


BOUNDARY CODE
  Inlet holds its value (by design). Surface reflection present.
