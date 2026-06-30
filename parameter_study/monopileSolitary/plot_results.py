"""
waves2Foam monopile analysis plots
Case: monopileSolitary
Run from the case root directory:
    python3 plot_results.py
Outputs: wave_gauges.png, cylinder_force.png, transmission.png
"""

import numpy as np
import matplotlib.pyplot as plt
import os

GAUGE_FILE  = "postProcessing/waveGauges/0/alpha.water"
FORCE_FILE  = "postProcessing/forces/0/force.dat"
CASE_NAME   = "monopileSolitary  (waves2Foam, solitaryFirst, H=0.10m, d=0.40m)"
GAUGE_X     = [1.0, 2.5, 4.0, 6.0, 7.5, 9.0]
STILL_WATER = 0.4
DOMAIN_H    = 0.7
CYLINDER_X  = 5.0
CYLINDER_R  = 0.5

def load_gauges(path):
    times, rows = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = line.split()
            times.append(float(vals[0]))
            rows.append([float(v) for v in vals[1:]])
    return np.array(times), np.array(rows).T

def alpha_to_eta(alpha):
    return STILL_WATER + (alpha - 0.5) * DOMAIN_H

def load_forces(path):
    times, fx, fy = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = line.split()
            times.append(float(vals[0]))
            fx.append(float(vals[1]))
            fy.append(float(vals[2]))
    return np.array(times), np.array(fx), np.array(fy)

def wave_height(eta, times, t_start=1.0, t_end=11.0):
    mask = (times >= t_start) & (times <= t_end)
    if mask.sum() == 0:
        return np.nan
    seg = eta[mask]
    return np.max(seg) - np.min(seg)

print("Loading gauge data ...")
t_g, alpha = load_gauges(GAUGE_FILE)
eta = alpha_to_eta(alpha)

colors = ['#1f4e79','#2e75b6','#9dc3e6','#833c00','#c55a11','#f4b183']
labels = [f"G{i+1}: x={GAUGE_X[i]}m ({'up' if i<3 else 'down'}stream)" for i in range(6)]

fig1, axes = plt.subplots(6, 1, figsize=(12, 14), sharex=True)
fig1.suptitle(f"Free Surface Elevation at Wave Gauges\n{CASE_NAME}", fontsize=11)
for i, ax in enumerate(axes):
    ax.plot(t_g, eta[i], color=colors[i], lw=1.2, label=labels[i])
    ax.axhline(STILL_WATER, color='gray', lw=0.7, ls='--')
    ax.set_ylabel("η (m)", fontsize=8)
    ax.legend(loc='upper right', fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(STILL_WATER - 0.08, STILL_WATER + 0.15)
axes[-1].set_xlabel("Time (s)", fontsize=10)
plt.tight_layout()
fig1.savefig("wave_gauges.png", dpi=150, bbox_inches='tight')
print("  -> wave_gauges.png")

print("Loading force data ...")
t_f, fx, fy = load_forces(FORCE_FILE)
peak_idx = np.argmax(np.abs(fx))
peak_t, peak_fx = t_f[peak_idx], fx[peak_idx]

fig2, ax2 = plt.subplots(figsize=(12, 4))
ax2.plot(t_f, fx, color='#c00000', lw=1.4, label='Fx (inline)')
ax2.plot(t_f, fy, color='#7030a0', lw=0.8, ls='--', label='Fy (transverse)', alpha=0.7)
ax2.axhline(0, color='black', lw=0.6)
ax2.axvline(peak_t, color='red', lw=0.8, ls=':', alpha=0.8)
ax2.annotate(f"Peak: {peak_fx:.2f} N\nt = {peak_t:.2f} s",
             xy=(peak_t, peak_fx), xytext=(peak_t+0.3, peak_fx*0.85),
             fontsize=8, color='#c00000',
             arrowprops=dict(arrowstyle='->', color='#c00000', lw=0.8))
ax2.set_xlabel("Time (s)", fontsize=10)
ax2.set_ylabel("Force (N)", fontsize=10)
ax2.set_title(f"Hydrodynamic Force on Cylinder\n{CASE_NAME}", fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
fig2.savefig("cylinder_force.png", dpi=150, bbox_inches='tight')
print("  -> cylinder_force.png")

H_gauges = [wave_height(eta[i], t_g) for i in range(6)]
H_inc  = np.nanmean(H_gauges[:3])
H_trans = np.nanmean(H_gauges[3:])
Kt = H_trans / H_inc if H_inc > 0 else np.nan

fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(13, 5))
fig3.suptitle(f"Wave Attenuation & Transmission\n{CASE_NAME}", fontsize=11)
ax3a.scatter(GAUGE_X[:3], H_gauges[:3], color='#2e75b6', s=80, zorder=5, label='Upstream')
ax3a.scatter(GAUGE_X[3:], H_gauges[3:], color='#c55a11', s=80, zorder=5, label='Downstream')
ax3a.axvspan(CYLINDER_X-CYLINDER_R, CYLINDER_X+CYLINDER_R, alpha=0.15, color='gray', label='Cylinder')
ax3a.set_xlabel("x (m)", fontsize=10)
ax3a.set_ylabel("Wave height H (m)", fontsize=10)
ax3a.set_title("Wave Height Along Flume", fontsize=10)
ax3a.legend(fontsize=8)
ax3a.grid(True, alpha=0.3)
bars = ax3b.bar(['Incident','Transmitted'], [H_inc, H_trans],
                color=['#2e75b6','#c55a11'], width=0.4, edgecolor='black', lw=0.8)
for bar, h in zip(bars, [H_inc, H_trans]):
    ax3b.text(bar.get_x()+bar.get_width()/2, h+0.001, f'{h:.4f} m',
              ha='center', va='bottom', fontsize=9)
ax3b.set_ylabel("Wave Height (m)", fontsize=10)
ax3b.set_title(f"Transmission Coefficient Kt = {Kt:.3f}", fontsize=10)
ax3b.set_ylim(0, max(H_inc, H_trans)*1.3)
ax3b.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
fig3.savefig("transmission.png", dpi=150, bbox_inches='tight')
print("  -> transmission.png")

print()
print("="*50)
print("RESULTS SUMMARY")
print("="*50)
for i in range(6):
    side = "upstream  " if i < 3 else "downstream"
    print(f"  G{i+1} x={GAUGE_X[i]:4.1f}m ({side}): H = {H_gauges[i]:.4f} m")
print(f"\nMean H upstream   : {H_inc:.4f} m")
print(f"Mean H downstream : {H_trans:.4f} m")
print(f"Transmission Kt   : {Kt:.3f}  ({Kt*100:.1f}% passes cylinder)")
print(f"\nPeak inline force : {peak_fx:.3f} N  at t = {peak_t:.3f} s")
print(f"Min inline force  : {fx.min():.3f} N  at t = {t_f[np.argmin(fx)]:.3f} s")
print("="*50)
