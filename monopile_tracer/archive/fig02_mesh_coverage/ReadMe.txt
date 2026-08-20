fig02_mesh_coverage
Where the LES mesh stops resolving the wake
===========================================

WHAT THIS SHOWS

  How much of the 600 x 120 m domain the LES mesh actually resolves.
  Binning the mesh onto a regular 0.5 m grid and counting how many cells
  received at least one mesh point gives a direct map of where the
  simulation has detail and where it does not.

  This is the diagnostic behind the caveat on every other figure here.

RESULT

  Coverage holds near 100% out to about x = 280 m, then collapses:

      x   0-100      90.6% of cells contain mesh points
      x 100-200      90.8%
      x 200-300      85.8%
      x 300-400      15.0%
      x 400-500       5.5%
      x 500-600       3.6%

  Across the full domain only 48.6% of cells are filled.

  Every filled column has all 84 of its vertical levels populated, none
  partial. The mesh is vertically extruded, so a horizontal position
  either has mesh points at every depth or none at all. That means finer
  vertical binning cannot recover the empty regions.

  In y, coverage is 53% through the middle band and drops to 44% at both
  edges, so the refined region spans a fixed lateral band rather than
  following the wake.

  Practical consequence: past about x = 290 m the mesh spacing is 1 to
  3 m. No binning recovers detail there, because the LES did not compute
  it. Coarser bins only average the same sparse data over a wider area.

INPUTS

  data/velocity_full_domain_05m/  -> symlink to shared data, 128 MB
      A single binned timestep at 0.5 m over the full domain, produced by
      monopile_tracer/scripts/extract_full.py. Only one timestep is
      needed, since the mesh is static.

TO REPRODUCE

  From this directory:

      python scripts/mesh_coverage.py    coverage by x and y band
      python scripts/mesh_bounds.py      finds where coverage drops off
      python scripts/plot_coverage.py    writes the coverage map

  Move the output PNG into results/ afterwards. Each runs in seconds.

NOTE

  The velocity file used is vel_595, from the earlier 399 to 595 s block.
  That does not matter here, because mesh coverage is a property of the
  mesh and does not change with time.
