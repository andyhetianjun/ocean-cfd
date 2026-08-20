task3_current_nocylinder
Regular wave on a uniform current, no structure
===============================================

PURPOSE

  Adds a uniform current to the task 1 setup, with no cylinder. This is the
  reference case for task 4: subtracting it isolates what the structure does
  when both waves and a current are present.

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

  Mesh            200 x 40 x 75 = 600,000 cells, identical to task 1.
                  Three vertical blocks refined through the free surface:
                      0.00 to 0.65 m    20 cells,  0.0325 m
                      0.65 to 0.95 m    40 cells,  0.0075 m
                      0.95 to 1.20 m    15 cells,  0.0167 m

  Boundaries      inlet         gabcVelocity, wave generation
                  outlet        gabcVelocity, absorption with current
                  bottom        slip
                  frontAndBack  slip
                  atmosphere    pressureInletOutletVelocity

  Run             0 to 12 s, adjustable timestep, Courant limit 0.25,
                  output every 0.25 s, Tsoft = 2 s.

  Gauge           single, x = 10 m, y = 2 m (flume centre).
                  One gauge is enough here, since with no cylinder there is
                  nothing to be upstream or downstream of.

RESULTS

  task3_surface_elevation_single_gauge.png    elevation at x = 10 m
  task3_surface_elevation_animation.gif       2-D surface elevation
  task3_velocity_animation.gif                2-D surface current velocity

TO REPRODUCE

  1. Mesh, initialise, run:
         blockMesh
         setWaveField
         decomposePar
         mpirun -np <N> waveFoam -parallel
         reconstructPar
  2. Gauge plot:
         python plot_elevation.py
  3. Animations: export PNG frames from ParaView into surface/ and
     velocity/, then
         convert -delay 4 -loop 0 surface/frames.*.png task3_surface_elevation_animation.gif

     Following make_gifMovie_T.sh from Yongxing Ma. -delay 4 gives 25 fps.

NOTES

  The current is following, in the same direction as the wave, which
  lengthens the apparent wavelength and lowers the encounter frequency
  relative to task 1.

  Slip on the bottom means the current has no boundary layer and is uniform
  with depth. Physically that is a simplification, but it matches the
  intent of a uniform current.
