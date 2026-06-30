#!/usr/bin/env python3
"""
Plot surface elevation at the single x=10m gauge with SWL and +/-A reference lines.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = "postProcessing/surfaceElevation/0/surfaceElevation.dat"
OUT_FILE = "task1_surface_elevation_single_gauge.png"

SWL = 0.80      # still water level (domain depth d)
H = 0.05        # wave height
A = H / 2.0     # amplitude

# Load data, skipping the header line and the 3 metadata lines (gauge x,y,z)
time = []
eta = []
with open(DATA_FILE) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("Time"):
            continue
        parts = line.split()
        t = float(parts[0])
        if t < 0:
            continue  # skip metadata rows (-1, -2, -3)
        time.append(t)
        eta.append(float(parts[1]))

time = np.array(time)
eta = np.array(eta)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(time, eta, color="tab:blue", linewidth=1.2, label="Surface elevation (x=10m)")
ax.axhline(SWL, color="black", linewidth=0.8, linestyle="-", label="SWL")
ax.axhline(SWL + A, color="tab:red", linewidth=0.8, linestyle="--", label=f"+A ({SWL+A:.3f} m)")
ax.axhline(SWL - A, color="tab:red", linewidth=0.8, linestyle="--", label=f"-A ({SWL-A:.3f} m)")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Surface elevation (m)")
ax.set_title("Task 1: Regular Wave Surface Elevation at x=10m (gauge_center)")
ax.legend(loc="upper right", fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, time.max())

plt.tight_layout()
plt.savefig(OUT_FILE, dpi=150)
print(f"Saved {OUT_FILE}")
print(f"Min eta: {eta.min():.4f}, Max eta: {eta.max():.4f}")
