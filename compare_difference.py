#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import os

TASK1_DIR = "/shared_folder/andyhe/project/waves2foam/task1_regularFlume/postProcessing/sampleSurface"
TASK2_DIR = "/shared_folder/andyhe/project/waves2foam/task2_cylinder_H005/postProcessing/sampleSurface"
OUT_DIR = "/shared_folder/andyhe/project/waves2foam/diff_frames"
os.makedirs(OUT_DIR, exist_ok=True)

xi = np.linspace(0, 20, 200)
yi = np.linspace(0, 4, 40)
XX, YY = np.meshgrid(xi, yi)

cyl_x, cyl_y, cyl_r = 5.0, 2.0, 0.25
SWL = 0.80

VMIN, VMAX = -0.025, 0.025
LEVELS = np.linspace(VMIN, VMAX, 51)

times = sorted([d for d in os.listdir(TASK1_DIR) if d.replace('.','').isdigit()], key=float)
print(f"Found {len(times)} timesteps")

for t in times:
    f1 = os.path.join(TASK1_DIR, t, "alpha.water_freeSurface.raw")
    f2 = os.path.join(TASK2_DIR, t, "alpha.water_freeSurface.raw")
    if not os.path.exists(f1) or not os.path.exists(f2):
        print(f"Skipping t={t}")
        continue

    d1 = np.loadtxt(f1, comments="#")
    d2 = np.loadtxt(f2, comments="#")

    eta1 = griddata((d1[:,0], d1[:,1]), d1[:,2] - SWL, (XX, YY), method='linear')
    eta2 = griddata((d2[:,0], d2[:,1]), d2[:,2] - SWL, (XX, YY), method='linear')

    diff = eta2 - eta1

    fig, ax = plt.subplots(figsize=(12, 3))
    cf = ax.contourf(XX, YY, diff, levels=LEVELS, cmap='RdBu_r', vmin=VMIN, vmax=VMAX, extend='both')
    cbar = plt.colorbar(cf, ax=ax, ticks=[-0.025, -0.0125, 0, 0.0125, 0.025])
    cbar.set_label('Δ Surface Elevation (m)')

    theta = np.linspace(0, 2*np.pi, 100)
    ax.fill(cyl_x + cyl_r*np.cos(theta), cyl_y + cyl_r*np.sin(theta), 'k')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(f'Wave Flume  t = {float(t):.2f} s  —  Difference (cylinder − no cylinder)')
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 4)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f"diff_{float(t):06.2f}.png"), dpi=100)
    plt.close()
    print(f"Saved t={t}")

print("Done!")
