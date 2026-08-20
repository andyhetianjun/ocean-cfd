fig03_domain_velocity
Cross-stream velocity over the full domain
==========================================

WHAT THIS SHOWS

  A single snapshot of depth-averaged cross-stream velocity across the
  whole 600 x 120 m domain. The alternating red and blue lobes downstream
  of the monopile are shed vortices, so this is the clearest single image
  of the vortex street.

  A dashed box marks the region the tracer animations cover.

RESULT

  results/full_domain_velocity.png

  The organised vortex street is visible from the pile out to roughly
  x = 300 m, then breaks down into faint streaks. That breakdown is the
  mesh coarsening, not the wake ending. See fig02.

  Two checks the script prints:

      empty columns      0.5145   fraction of the grid with no mesh points
      upstream u         0.500    against an expected 0.496 m/s

  The second confirms the depth-averaging and the counts weighting are
  correct, since undisturbed flow upstream of the pile should sit at the
  inlet value.

INPUTS

  data/velocity_full_domain_05m/  -> symlink to shared data, 128 MB
      Single binned timestep at 0.5 m over the full domain, from
      monopile_tracer/scripts/extract_full.py.

TO REPRODUCE

  From this directory:

      python scripts/plot_full_domain.py

  Runs in seconds. Writes the PNG to the working directory, so move it
  into results/ afterwards.

NOTE

  This snapshot is timestep 595, from the earlier 399 to 595 s block,
  while the tracer animations in fig01 use 1049 to 1400 s. The vortex
  street looks the same either way, but the two are not the same instant.
  To make the set consistent, point the script at a timestep in
  monopile_tracer/data/velocity_full instead.
