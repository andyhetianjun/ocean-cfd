FULL-DOMAIN SLICK, WITH BUOYANCY
================================
See ../../ReadMe.txt for the velocity source, boundary conditions, cache
mechanics, and the meaning of the two limits.

Oil covers the entire 600 x 120 m surface at t = 0, with the inlet
holding that value so the sheet is continuous rather than draining.

PARAMETERS
  OIL_MASS_KG      573.0    over 72000 m2
  OIL_FULL_DOMAIN  True
  VERTICAL_MIX     False    w zeroed at load
  TRACER_NZ        100      dz = 0.399 m

  OIL_RADIUS and OIL_OFFSET are unused when OIL_FULL_DOMAIN is set. So is
  BEARING, though it still appears in the log header - ignore it there.

DERIVED
  Areal loading   7.96 g/m2   (chosen to be comparable to the disc cases)
  Film thickness  9.4 um
  Surface cell    19.9 g/m3

RESULT
  results/oil_fulldomain_buoyant.gif
  Colour range 7.958 to 7.966 g/m2.

This case shows no signal, and that is the correct result rather than a
failure. With w = 0 the oil moves only horizontally, and translating a
horizontally uniform sheet leaves it uniform - there is no edge to
distort and no gradient to stretch.

The colour bar therefore spans 0.008 g/m2, a zoom into numerical noise at
the 0.1 percent level. The faint band near the inlet is numerical
diffusion from the boundary condition propagating downstream at the flow
speed; it reaches about 176 m by the end of the run, which matches 352 s
at 0.5 m/s. It is not physical structure.

Worth stating explicitly whenever this figure is shown, or it looks like
something went wrong.

Mass holds flat at 573 kg to within 0.1 percent - inflow and outflow
balance for a uniform sheet in a uniform mean flow.

Note TRACER_NZ is 100 here against 150 in ../w_on. With w = 0 the
vertical grid does nothing, so the coarser grid saves time and changes
nothing. It does mean the pair is not strictly like for like.

BOUNDARY CODE
  Inlet holds its value (by design). No surface reflection - nothing to
  reflect with w = 0.
