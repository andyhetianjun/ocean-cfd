waves2Foam wave flume cases
===========================

WHAT THIS IS

  Four OpenFOAM cases run with the waves2Foam toolbox, built to learn the
  toolbox and to document a working setup so it can be picked up again.

  The four form a 2 x 2. Waves are present in all of them; the current and
  the cylinder are switched on and off:

                        no cylinder            cylinder
      no current        task1_regularFlume     task2_cylinder
      current 0.2 m/s   task3_current_...      task4_cylinder_current

  That means each effect can be isolated by subtracting one case from
  another. task1 vs task2 gives the cylinder's effect on waves alone;
  task3 vs task4 gives it with a current present.

SHARED SETUP

  Every case uses the same flume and the same wave. Only the current and
  the cylinder change.

      Solver        waveFoam (waves2Foam on OpenFOAM v2412), laminar
      Flume         20 m long, 4 m wide, 1.2 m tall
      Water depth   0.8 m
      Wave          Stokes 1st order, H = 0.05 m, T = 3 s, L = 7.90 m
      Mesh          200 x 40 x 75 = 600,000 cells, refined vertically
                    through the free surface
      Boundaries    GABC at inlet and outlet, slip on bottom and sides
      Run           12 s (four wave periods), Courant limit 0.25

  Full details, including gauge positions and cylinder geometry, are in
  each case's own ReadMe.txt.

CONTENTS

  task1_regularFlume/         waves only, empty flume
  task2_cylinder/             waves + cylinder
  task3_current_nocylinder/   waves + current
  task4_cylinder_current/     waves + current + cylinder

  compare_difference.py       builds task1 vs task2 difference animation
  compare_difference_t3_t4.py builds task3 vs task4 difference animation
  downloadFromNibi.sh         pulls case output from the Nibi cluster

  parameter_study/            wave parameter sweep, see sweep_results.txt
  legacy_olaflow_monopile/    earlier olaFlow attempt, kept for reference
  regularWave/, squarePile/   reference cases, not tracked

WORKFLOW

  Each case follows the same path:

      blockMesh                 build the base mesh
      snappyHexMesh -overwrite  cut in the cylinder (tasks 2 and 4 only)
      setWaveField              initialise the wave field
      decomposePar              split for parallel running
      mpirun -np <N> waveFoam -parallel
      reconstructPar            reassemble

  Surface elevation is recorded by the surfaceElevation function object in
  system/controlDict, with gauge positions in waveGaugesNProbes/. The
  plot_elevation*.py scripts turn that output into time series plots.

  Animations are exported as PNG frame sequences from ParaView and
  assembled with ImageMagick:

      convert -delay 4 -loop 0 velocity/frames.*.png output.gif

  following make_gifMovie_T.sh from Yongxing Ma.

NOT TRACKED

  Time directories, processor* folders, generated mesh, postProcessing/ and
  ParaView frame exports are all regenerable and excluded via .gitignore.
  The finished animations and plots are tracked.

WORTH KNOWING

  Absorption uses generalised absorbing boundary conditions rather than
  relaxation zones, so relaxationNames is empty in waveProperties. The GABC
  polynomial is fitted for kh from 0 to 3.0 and these cases run at kh = 0.64.

  Bottom and sidewalls are slip, so there are no boundary layers anywhere.
  Fine for wave propagation, but these cases say nothing about bed shear,
  and the current is uniform with depth.

  Blockage is D/W = 12.5%, high enough that some sidewall interaction with
  the cylinder is plausible.

  The ParaView pipeline that produced the animations was built interactively
  and has not been saved. Recording it with Tools > Start Trace, or saving a
  .pvsm state file, would close the last gap in reproducing these results.
