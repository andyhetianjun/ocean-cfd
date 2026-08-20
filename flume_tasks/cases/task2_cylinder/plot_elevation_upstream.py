#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = "postProcessing/surfaceElevation/0/surfaceElevation.dat"
OUT_FILE = "task2_H005_surface_elevation_upstream.png"

SWL = 0.80
H = 0.05
A = H / 2.0

time = []
eta_up = []

with open(DATA_FILE) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("Time"):
            continue
        parts = line.split()
        t = float(parts[0])
        if t < 0:
            continue
        time.append(t)
        eta_up.append(float(parts[1]))

time = np.array(time)
eta_up = np.array(eta_up) - SWL

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(time, eta_up, color="tab:blue", linewidth=1.2, label="Upstream (x=3m)")
ax.axhline( A, color="tab:red", linewidth=0.8, linestyle="--", label=f"+A ({A:.3f} m)")
ax.axhline(-A, color="tab:red", linewidth=0.8, linestyle="--", label=f"-A ({-A:.3f} m)")
ax.axhline(0, color="black", linewidth=0.6, linestyle="-")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Surface Elevation (m)")
ax.set_title("Task 2 (H=0.05m): Surface Elevation Upstream of Cylinder (x=3m)")
ax.legend(loc="upper right", fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, time.max())
plt.tight_layout()
plt.savefig(OUT_FILE, dpi=150)
print(f"Saved {OUT_FILE}")
