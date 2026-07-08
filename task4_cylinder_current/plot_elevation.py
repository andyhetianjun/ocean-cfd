#!/usr/bin/env python3
"""
Plot surface elevation deviation from SWL at upstream and downstream gauges.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = "postProcessing/surfaceElevation/0/surfaceElevation.dat"
OUT_FILE = "task2_surface_elevation.png"

SWL = 0.80
H = 0.10
A = H / 2.0

time = []
eta_up = []
eta_down = []

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
        eta_down.append(float(parts[2]))

time = np.array(time)
eta_up = np.array(eta_up) - SWL
eta_down = np.array(eta_down) - SWL

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(time, eta_up, color="tab:blue", linewidth=1.2, label="Upstream (x=3m)")
ax.plot(time, eta_down, color="tab:orange", linewidth=1.2, label="Downstream (x=7m)")
ax.axhline( A, color="tab:red", linewidth=0.8, linestyle="--", label=f"+A ({A:.3f} m)")
ax.axhline(-A, color="tab:red", linewidth=0.8, linestyle="--", label=f"-A ({-A:.3f} m)")
ax.axhline(0, color="black", linewidth=0.6, linestyle="-")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Surface elevation deviation from SWL (m)")
ax.set_title("Task 2: Regular Wave Surface Elevation at Upstream and Downstream Gauges")
ax.legend(loc="upper right", fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, time.max())

plt.tight_layout()
plt.savefig(OUT_FILE, dpi=150)
print(f"Saved {OUT_FILE}")
print(f"Upstream   - Min: {eta_up.min():.4f}, Max: {eta_up.max():.4f}")
print(f"Downstream - Min: {eta_down.min():.4f}, Max: {eta_down.max():.4f}")
