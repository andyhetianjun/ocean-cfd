task4_cylinder_current
Regular wave on a uniform current, passing a cylinder
=====================================================

PURPOSE

  Waves, a uniform current, and a surface-piercing cylinder together. This
  was the case the task list was originally building toward; task 3 was
  added afterwards as the no-structure reference so the two could be
  subtracted.

MODEL DESCRIPTION

  Solver          waveFoam (waves2Foam on OpenFOAM v2412), laminar

  Flume           20 m long, 4 m wide, 1.2 m tall
  Water depth     0.8 m
  Wave theory     Stokes 1st order
  Wave height     0.05 m
  Wave period     3 s
  Wavelength      7.90 m  (k = 0.795 1/m, omega = 2.094 rad/s, kh = 0.64)
  Current         0.2 m/s in +x, entering at the left boundary, set through
                  the outlet potentialCurrent U = (0.2 0 0)

  Cylinder        diameter 0.5 m, full height 0 to 1.2 m, surface piercing
                  Centred at x = 5 m, y = 2 m
                  Geometry in constant/triSurface/cylinder.stl

                  D/L = 0.063 and kD = 0.40, so inertia dominated with weak
                  diffraction. Blockage D/W = 12.5%.

  Mesh            Base 200 x 40 x 75 = 600,000 cells.
                  snappyHexMesh refines the cylinder surface to level 2-3,
                  roughly 0.0125 m at the pile, with a refinement box from
                  (4.0, 1.25, 0.0) to (6.0, 2.75, 1.2) and
                  nCellsBetweenLevels 3.

  Boundaries      inlet         gabcVelocity, wave generation
                  outlet        gabcVelocity, absorption with current
                  bottom        slip
                  frontAndBack  slip
                  atmosphere    pressureInletOutletVelocity
                  cylinder      noSlip

  Run             0 to 12 s, adjustable timestep, Courant limit 0.25,
                  output every 0.25 s, Tsoft = 2 s.

  Gauges          upstream    x = 3 m, y = 2 m
                  downstream  x = 7 m, y = 2 m
                  Column 1 is upstream, column 2 downstream in
                  surfaceElevation.dat.

RESULTS

  task4_surface_elevation_upstream.png     elevation at x = 3 m
  task4_surface_elevation_downstream.png   elevation at x = 7 m
  task4_surface_elevation_animation.gif    2-D surface elevation
  task4_velocity_animation.gif             2-D surface current velocity

  ../task3_vs_task4_difference.gif isolates the cylinder's effect by
  subtracting task 3, built by ../compare_difference_t3_t4.py.

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
         convert -delay 4 -loop 0 surface/frames.*.png task4_surface_elevation_animation.gif

     Following make_gifMovie_T.sh from Yongxing Ma. -delay 4 gives 25 fps.

NOTES

  This case differs from task 2 only by the 0.2 m/s current, and from task 3
  only by the cylinder. That makes both single-variable comparisons
  available.

  The current is following, so wave and current interact rather than oppose.
  Combined with the no-slip cylinder in a slip-walled flume, the wake here
  is driven by the structure alone with no bed or wall boundary layers.
