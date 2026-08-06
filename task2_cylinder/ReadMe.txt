task2_cylinder
Regular wave passing a surface-piercing cylinder
================================================

PURPOSE

  Same wave and flume as task 1, with a vertical cylinder added. Comparing
  the two isolates what the structure does to the wave field.

MODEL DESCRIPTION

  Solver          waveFoam (waves2Foam on OpenFOAM v2412), laminar

  Flume           20 m long, 4 m wide, 1.2 m tall
  Water depth     0.8 m
  Wave theory     Stokes 1st order
  Wave height     0.05 m
  Wave period     3 s
  Wavelength      7.90 m  (k = 0.795 1/m, omega = 2.094 rad/s, kh = 0.64)
  Current         none

  Cylinder        diameter 0.5 m, full height 0 to 1.2 m, so surface piercing
                  Centred at x = 5 m, y = 2 m (flume centre laterally,
                  5 m from the inlet, 0.63 wavelengths downstream)
                  Geometry in constant/triSurface/cylinder.stl

                  D/L  = 0.063 and kD = 0.40, so the cylinder is small
                  relative to the wave. Inertia dominated, weak diffraction.
                  Blockage D/W = 12.5%, high enough that some sidewall
                  interaction is plausible.

  Mesh            Base 200 x 40 x 75 = 600,000 cells, as task 1.
                  snappyHexMesh refines the cylinder surface to level 2-3,
                  giving roughly 0.0125 m at the pile, with a refinement box
                  from (4.0, 1.25, 0.0) to (6.0, 2.75, 1.2) and
                  nCellsBetweenLevels 3.

  Boundaries      inlet         gabcVelocity, wave generation
                  outlet        gabcVelocity, wave absorption
                  bottom        slip
                  frontAndBack  slip
                  atmosphere    pressureInletOutletVelocity
                  cylinder      noSlip

  Run             0 to 12 s, adjustable timestep, Courant limit 0.25,
                  output every 0.25 s, Tsoft = 2 s.

  Gauges          upstream    x = 3 m, y = 2 m  (2 m before the pile)
                  downstream  x = 7 m, y = 2 m  (2 m after the pile)
                  Both vertical lines sampled at 100 points, 0 to 1.2 m.
                  In surfaceElevation.dat, column 1 is upstream and
                  column 2 is downstream.

RESULTS

  task2_H005_surface_elevation_upstream.png     elevation at x = 3 m
  task2_H005_surface_elevation_downstream.png   elevation at x = 7 m
  task2_velocity_animation.gif                  2-D surface current velocity
  task2_H005_velocity_animation.gif             velocity, H = 0.05 m run
  task2_H005_birdseye_elevation.gif             surface elevation, plan view
  task2_difference_animation.gif                task 2 minus task 1

  ../task1_vs_task2_difference.gif is the same comparison at the repository
  root, built by ../compare_difference.py.

TO REPRODUCE

  1. Mesh, initialise, run:
         blockMesh
         snappyHexMesh -overwrite
         setWaveField
         decomposePar
         mpirun -np <N> waveFoam -parallel
         reconstructPar
  2. Gauge plots:
         python plot_elevation_upstream.py
         python plot_elevation_downstream.py
  3. Animations: export PNG frames from ParaView, then
         convert -delay 4 -loop 0 velocity/frames.*.png task2_velocity_animation.gif

     Following make_gifMovie_T.sh from Yongxing Ma. -delay 4 gives 25 fps.

NOTES

  The H005 prefix on some outputs refers to the 0.05 m wave height, from
  when several heights were being compared. All results here are H = 0.05 m.

  Both gauges sit on the centreline, so the downstream one is directly in
  the cylinder's shadow.
