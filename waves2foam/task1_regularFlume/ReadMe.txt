task1_regularFlume
Regular wave propagation in an empty wave flume
===============================================

PURPOSE

  Baseline. A regular wave travelling down an empty flume, no structure and
  no current. Everything else in this set is measured against it, and it was
  the first case built while learning waves2Foam.

MODEL DESCRIPTION

  Solver          waveFoam (waves2Foam on OpenFOAM v2412), laminar

  Flume           20 m long, 4 m wide, 1.2 m tall
  Water depth     0.8 m, leaving 0.4 m of air above
  Wave theory     Stokes 1st order
  Wave height     0.05 m
  Wave period     3 s
  Wavelength      7.90 m  (k = 0.795 1/m, omega = 2.094 rad/s, kh = 0.64)
  Current         none

  Mesh            200 x 40 x 75 = 600,000 cells, uniform in x and y (0.1 m)
                  Three vertical blocks, refined through the free surface:
                      0.00 to 0.65 m    20 cells,  0.0325 m
                      0.65 to 0.95 m    40 cells,  0.0075 m
                      0.95 to 1.20 m    15 cells,  0.0167 m
                  About 79 cells per wavelength horizontally.

  Boundaries      inlet         gabcVelocity, wave generation
                  outlet        gabcVelocity, wave absorption
                  bottom        slip, no bottom friction
                  frontAndBack  slip, no sidewall friction
                  atmosphere    pressureInletOutletVelocity

                  Absorption uses generalised absorbing boundary conditions
                  rather than relaxation zones, which is why relaxationNames
                  is empty. The GABC polynomial is fitted for kh from 0 to
                  3.0; this case runs at kh = 0.64, well inside that.

  Run             0 to 12 s (four wave periods), adjustable timestep,
                  Courant limit 0.25, output every 0.25 s.
                  Tsoft = 2 s ramps the wave in to avoid a startup shock.

  Gauge           single, x = 10 m, y = 2 m (flume centre)

RESULTS

  task1_surface_elevation_single_gauge.png   elevation time series
  task1_velocity_animation.gif               2-D surface current velocity
  task1_birdseye_animation_v4.gif            2-D surface elevation, plan view

TO REPRODUCE

  1. Mesh and initialise:
         blockMesh
         setWaveField
  2. Run:
         decomposePar
         mpirun -np <N> waveFoam -parallel
         reconstructPar
  3. Gauge plot:
         python plot_elevation.py
  4. Animations: export PNG frames from ParaView into surface/ and
     velocity/, then assemble with ImageMagick:
         convert -delay 4 -loop 0 velocity/frames.*.png task1_velocity_animation.gif

     The convert command follows make_gifMovie_T.sh from Yongxing Ma.
     -delay 4 gives 25 fps.

  Time directories, processor* folders and postProcessing/ are not tracked,
  since they regenerate from the above.

NOTES

  Twelve seconds is four wave periods. Enough to see steady propagation once
  the 2 s ramp has passed, not enough for long-term statistics.

  Slip on the bottom and sidewalls means no boundary layers anywhere. Fine
  for wave propagation, but this case says nothing about bed shear.
