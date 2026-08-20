VORTICITY-Z ANIMATION
=====================

Two-panel animation of the vertical component of vorticity in the wake of
the monopile.

  omega_z = dv/dx - du/dy      units 1/s

  top panel     plan view at mid-depth, x = 0-300 m, y = 0-120 m
  bottom panel  vertical section along the pile centreline (y = 60 m),
                depth 0 at the surface down to about -40 m at the bed

DATA SOURCE
-----------
Reads Andy's binned extraction of Yongxing's LES:

  /shared_folder/andyhe/project/monopile_tracer/data/velocity_full_025

Regular 0.25 m grid, 2400 x 480 x 84, timesteps 1049-1400 s (352 steps),
produced from the raw OpenFOAM output at

  /shared_folder/yongxing/OpenFOAM/simulationCases/flow_past_cylinder/
      realistic_cases/domain_D8X600Y120/uniform

Domain 600 x 120 x 40 m, monopile D = 8 m at (40, 60) spanning the full
depth, freestream about 0.5 m/s, slip bed.

Two limitations worth knowing:

- The extraction fills roughly 60% of cells by nearest neighbour.
  Vorticity is a spatial derivative, so it amplifies that. The result is
  grainier than one computed on the native mesh.

- The LES mesh is refined only to about x = 290 m and coarsens to 1-3 m
  beyond. The figure is cropped at x = 300 for that reason.

RUNNING
-------
Run from the folder containing data/, NOT from inside scripts/ - paths
are relative.

  cd vorticity
  /home/andyhe/tracer_env/bin/python -u scripts/vorticity_z_v2.py

Or detached:

  setsid nohup /home/andyhe/tracer_env/bin/python -u \
      scripts/vorticity_z_v2.py > render.log 2>&1 < /dev/null &

  tail -f render.log

The shell prints "Done" immediately with setsid - that means detached,
not finished. Check with:

  pgrep -af vorticity_z_v2.py

Output: results/vorticity_z_2panel.gif

About 35 minutes for 352 frames. Memory stays under 3 GB - it loads one
timestep, takes the two slices, discards the rest.

TESTING FIRST
-------------
  cp scripts/vorticity_z_v2.py scripts/test.py
  sed -i 's|range(1049, 1401)|range(1049, 1059)|' scripts/test.py
  /home/andyhe/tracer_env/bin/python -u scripts/test.py

10 frames, under a minute.

SETTINGS
--------
All near the top of the script.

  VLIM = 0.15     Colour limits, symmetric, 1/s. The 99th percentile of
                  |omega_z| at mid-depth is about 0.07. Drop to 0.08 if
                  the wake looks pale, raise if saturated.

  X_MAX = 300.0   Downstream crop. Don't push past 300 - mesh coarsens.

  PLAN_Z = 0.0    Plan view depth in model coords, where the surface is
                  +19.76 and the bed -19.76. So 0.0 is mid-depth. Use
                  19.0 for near-surface.

  SECT_Y = 60.0   Section location. 60 is the pile centreline.

  FPS = 24        GIF stores delay in hundredths of a second, so 24 fps
                  is actually written as 25.

  height_ratios   In the GridSpec call, currently [2.0, 1.0]. Lower the
                  first number for a taller bottom panel.

NOTES ON THE FIGURE
-------------------
Red and blue are opposite signs of rotation. The alternating pattern
downstream is vortex shedding.

The depth axis runs 0 at the surface to -40 at the bed. Data actually
stops at -39.5, since the grid holds cell centres rather than faces, so
the -40 tick sits just outside the data.

Plan view is true 1:1 aspect. Cropped to 300 m gives a 2.5:1 panel.

SMALLER FILES
-------------
GIF is large - the earlier three-panel version was 51 MB, too big to
email. If ffmpeg is available, mp4 is far smaller, but that needs frames
written as PNGs first rather than straight to GIF.

SHARED MACHINE
--------------
Check free -g before starting anything large. The OOM killer picks the
largest process. This script is light, but the tracer runs in the
neighbouring folders need 20-30 GB.

Worth printing a running diagnostic every N frames. A bug in the tracer
code was caught purely because a printed mass total was climbing when it
should have been flat.

Questions to Andy.
