#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import os

TASK3_DIR = "/shared_folder/andyhe/project/waves2foam/task3_current_nocylinder/postProcessing/sampleSurface"
TASK4_DIR = "/shared_folder/andyhe/project/waves2foam/task4_cylinder_current/postProcessing/sampleSurface"
OUT_DIR = "/shared_folder/andyhe/project/waves2foam/diff_frames_t3_t4"
os.makedirs(OUT_DIR, exist_ok=True)

xi = np.linspace(0, 20, 200)
yi = np.linspace(0, 4, 40)
XX, YY = np.meshgrid(xi, yi)

cyl_x, cyl_y, cyl_r = 5.0, 2.0, 0.25
SWL = 0.80

VMIN, VMAX = -0.009, 0.009
LEVELS = np.linspace(VMIN, VMAX, 51)

times = sorted([d for d in os.listdir(TASK3_DIR) if d.replace('.','').isdigit()], key=float)
print(f"Found {len(times)} timesteps")

for t in times:
    f3 = os.path.join(TASK3_DIR, t, "alpha.water_freeSurface.raw")
    f4 = os.path.join(TASK4_DIR, t, "alpha.water_freeSurface.raw")
    if not os.path.exists(f3) or not os.path.exists(f4):
        print(f"Skipping t={t}")
        continue

    d3 = np.loadtxt(f3, comments="#")
    d4 = np.loadtxt(f4, comments="#")

    eta3 = griddata((d3[:,0], d3[:,1]), d3[:,2] - SWL, (XX, YY), method='linear')
    eta4 = griddata((d4[:,0], d4[:,1]), d4[:,2] - SWL, (XX, YY), method='linear')

    diff = eta4 - eta3

    fig, ax = plt.subplots(figsize=(12, 3))
    cf = ax.contourf(XX, YY, diff, levels=LEVELS, cmap='RdBu_r', vmin=VMIN, vmax=VMAX, extend='both')
    cbar = plt.colorbar(cf, ax=ax, ticks=[-0.009, -0.0045, 0, 0.0045, 0.009])
    cbar.set_label('Δ Surface Elevation (m)')

    theta = np.linspace(0, 2*np.pi, 100)
    ax.fill(cyl_x + cyl_r*np.cos(theta), cyl_y + cyl_r*np.sin(theta), 'k')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(f'Wave Flume  t = {float(t):.2f} s  —  Difference (Task 4 - Task 3)')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 4)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f"diff_{float(t):06.2f}.png"), dpi=100)
    plt.close()
    print(f"Saved t={t}")

print("Done!")
