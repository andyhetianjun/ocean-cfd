"""
Parameter sweep comparison — works from postProcessing/ only.
Run from: /shared_folder/andyhe/project/waves2foam/
python3 compare_cases.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

WORK    = "/shared_folder/andyhe/project/waves2foam"
SWL     = 0.4
DOMAIN_H = 0.7

CASES = [
    ("monopileSolitary", "solitaryFirst",  0.10, None, "Solitary\nH=0.10"),
    ("solitary_H005",    "solitaryFirst",  0.05, None, "Solitary\nH=0.05"),
    ("solitary_H015",    "solitaryFirst",  0.15, None, "Solitary\nH=0.15"),
    ("solitary_H020",    "solitaryFirst",  0.20, None, "Solitary\nH=0.20"),
    ("regular_T2",       "stokesFirst",    0.10, 2.0,  "Regular\nT=2s"),
    ("regular_T3",       "stokesFirst",    0.10, 3.0,  "Regular\nT=3s"),
    ("regular_T5",       "stokesFirst",    0.10, 5.0,  "Regular\nT=5s"),
    ("regular_T8",       "stokesFirst",    0.10, 8.0,  "Regular\nT=8s"),
    ("stokes2_T3",       "stokesSecond",   0.10, 3.0,  "Stokes2nd\nT=3s"),
    ("regular_T3_H005",   "stokesFirst",    0.05, 3.0,  "Regular\nH=0.05 T=3s"),
    ("regular_T3_H015",   "stokesFirst",    0.15, 3.0,  "Regular\nH=0.15 T=3s"),
    ("regular_T3_H020",   "stokesFirst",    0.20, 3.0,  "Regular\nH=0.20 T=3s"),
    ("chappelear_H005",   "chappelear1962", 0.05, None, "Chappelear\nH=0.05"),
    ("chappelear_H010",   "chappelear1962", 0.10, None, "Chappelear\nH=0.10"),
    ("chappelear_H015",   "chappelear1962", 0.15, None, "Chappelear\nH=0.15"),
    ("chappelear_H020",   "chappelear1962", 0.20, None, "Chappelear\nH=0.20"),
]

def load_alpha(case_dir):
    fpath = os.path.join(case_dir, "postProcessing/waveGauges/0/alpha.water")
    if not os.path.exists(fpath): return None, None
    times, rows = [], []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            vals = line.split()
            times.append(float(vals[0]))
            rows.append([float(v) for v in vals[1:]])
    return np.array(times), np.array(rows).T  # (6, n_times)

def load_forces(case_dir):
    fpath = os.path.join(case_dir, "postProcessing/forces/0/force.dat")
    if not os.path.exists(fpath): return None, None, None
    times, fx, fy = [], [], []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            vals = line.split()
            times.append(float(vals[0]))
            fx.append(float(vals[1]))
            fy.append(float(vals[2]))
    return np.array(times), np.array(fx), np.array(fy)


def load_surface_elevation(case_dir):
    fpath = os.path.join(case_dir, "postProcessing/surfaceElevation/0/surfaceElevation.dat")
    if not os.path.exists(fpath):
        return None, None
    times, rows = [], []
    with open(fpath) as f:
        for i, line in enumerate(f):
            if i < 4: continue
            line = line.strip()
            if not line: continue
            vals = line.split()
            times.append(float(vals[0]))
            rows.append([float(v) for v in vals[1:]])
    if not times:
        return None, None
    return np.array(times), np.array(rows).T

def alpha_to_eta(alpha):
    return SWL + (alpha - 0.5) * DOMAIN_H

def wave_height(eta, t, t0=1.5, t1=11.0):
    m = (t >= t0) & (t <= t1)
    if m.sum() < 5: return np.nan
    return eta[m].max() - eta[m].min()

# ── collect all results ───────────────────────────────────────────────────
results = []
print(f"{'Case':<22} {'H_inc':>7} {'H_trans':>8} {'Kt':>6} {'Fx_peak':>9} {'Fx_min':>9} {'t_peak':>7}")
print("-"*75)

for (name, wtype, H, T, label) in CASES:
    cdir = os.path.join(WORK, name)
    if not os.path.isdir(cdir):
        print(f"{name:<22}  -- not found --")
        continue

    t_a, alpha = load_alpha(cdir)
    t_f, fx, fy = load_forces(cdir)

    t_e, eta_e = load_surface_elevation(cdir)
    if t_e is not None:
        Hi = eta_e[1][4:].max() - 0.4  # G2 solitary peak above SWL
        Ht = ((eta_e[4][4:].max() - 0.4) + (eta_e[5][4:].max() - 0.4)) / 2
        Kt = Ht/Hi if (Hi and Hi > 0) else np.nan
    elif t_a is not None:
        eta = alpha_to_eta(alpha)
        Hi  = wave_height(eta[1], t_a)
        Ht  = np.nanmean([wave_height(eta[4], t_a), wave_height(eta[5], t_a)])
        Kt  = Ht/Hi if (Hi and Hi > 0) else np.nan
    else:
        Hi = Ht = Kt = np.nan

    if t_f is not None:
        peak_idx = np.argmax(np.abs(fx))
        Fx_peak  = fx[peak_idx]
        Fx_min   = fx.min()
        t_peak   = t_f[peak_idx]
    else:
        Fx_peak = Fx_min = t_peak = np.nan

    results.append(dict(name=name, label=label, wtype=wtype,
                        H=H, T=T, Hi=Hi, Ht=Ht, Kt=Kt,
                        Fx_peak=Fx_peak, Fx_min=Fx_min, t_peak=t_peak,
                        t_f=t_f, fx=fx, t_a=t_a, alpha=alpha if t_a is not None else None))

    Hi_s  = f"{Hi:.4f}"  if not np.isnan(Hi)     else "  nan "
    Ht_s  = f"{Ht:.4f}"  if not np.isnan(Ht)     else "   nan  "
    Kt_s  = f"{Kt:.3f}"  if not np.isnan(Kt)     else "  nan"
    Fp_s  = f"{Fx_peak:.2f}" if Fx_peak and not np.isnan(Fx_peak) else "   nan  "
    Fm_s  = f"{Fx_min:.2f}"  if Fx_min  and not np.isnan(Fx_min)  else "   nan  "
    tp_s  = f"{t_peak:.3f}"  if t_peak  and not np.isnan(t_peak)  else "  nan"
    print(f"{name:<22} {Hi_s:>7} {Ht_s:>8} {Kt_s:>6} {Fp_s:>9} {Fm_s:>9} {tp_s:>7}")

# ── Figure 1: Force time series — solitary wave height comparison ─────────
fig, axes = plt.subplots(4, 1, figsize=(12, 12), sharex=True)
fig.suptitle("Inline Force on Cylinder — Solitary Wave Height Study\n(waves2Foam, solitaryFirst, d=0.40m)", fontsize=11)
sol_cases = [r for r in results if r["wtype"] == "solitaryFirst"]
colors = ['#1f4e79','#2e75b6','#c55a11','#c00000']
for i, (r, ax, col) in enumerate(zip(sol_cases, axes, colors)):
    if r["t_f"] is None: continue
    ax.plot(r["t_f"], r["fx"], color=col, lw=1.2,
            label=f"H={r['H']}m  Fx_peak={r['Fx_peak']:.1f}N")
    ax.axhline(0, color='black', lw=0.5)
    ax.set_ylabel("Fx (N)", fontsize=8)
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
axes[-1].set_xlabel("Time (s)", fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(WORK, "fig1_solitary_force_timeseries.png"), dpi=150, bbox_inches='tight')
print("\nSaved: fig1_solitary_force_timeseries.png")

# ── Figure 2: Peak force vs wave height (solitary) ────────────────────────
sol = [r for r in results if r["wtype"]=="solitaryFirst" and r["Fx_peak"] and not np.isnan(r["Fx_peak"])]
sol_sorted = sorted(sol, key=lambda r: r["H"])
if sol_sorted:
    Hs  = [r["H"] for r in sol_sorted]
    Fps = [r["Fx_peak"] for r in sol_sorted]
    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot(Hs, Fps, 'o-', color='#c00000', lw=2, ms=9, zorder=5)
    for r in sol_sorted:
        ax.annotate(f"{r['Fx_peak']:.1f} N",
                    (r["H"], r["Fx_peak"]),
                    textcoords="offset points", xytext=(8,4), fontsize=9)
    # theoretical H^2 scaling reference
    H_ref = np.linspace(min(Hs), max(Hs), 50)
    F_ref = Fps[1] * (H_ref/Hs[1])**2
    ax.plot(H_ref, F_ref, '--', color='gray', lw=1, label='H² scaling (linear ref)')
    ax.set_xlabel("Wave Height H (m)", fontsize=11)
    ax.set_ylabel("Peak Inline Force Fx (N)", fontsize=11)
    ax.set_title("Peak Force vs Wave Height\n(solitaryFirst, d=0.40m, cylinder r=0.5m)", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(WORK, "fig2_force_vs_H.png"), dpi=150, bbox_inches='tight')
    print("Saved: fig2_force_vs_H.png")

# ── Figure 3: Force time series — regular wave period comparison ──────────
reg_cases = [r for r in results if r["wtype"] in ("stokesFirst","stokesSecond") and r["t_f"] is not None]
reg_cases = sorted(reg_cases, key=lambda r: r["T"] or 0)
if reg_cases:
    fig, axes = plt.subplots(len(reg_cases), 1, figsize=(12, 3*len(reg_cases)), sharex=True)
    if len(reg_cases)==1: axes=[axes]
    fig.suptitle("Inline Force on Cylinder — Regular Wave Study\n(H=0.10m, d=0.40m)", fontsize=11)
    colors_reg = ['#1f4e79','#2e75b6','#70ad47','#c55a11','#c00000']
    for i, (r, ax) in enumerate(zip(reg_cases, axes)):
        col = colors_reg[i % len(colors_reg)]
        ax.plot(r["t_f"], r["fx"], color=col, lw=1.0,
                label=f"{r['label'].replace(chr(10),' ')}  Fx_peak={r['Fx_peak']:.1f}N")
        ax.axhline(0, color='black', lw=0.5)
        ax.set_ylabel("Fx (N)", fontsize=8)
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("Time (s)", fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(WORK, "fig3_regular_force_timeseries.png"), dpi=150, bbox_inches='tight')
    print("Saved: fig3_regular_force_timeseries.png")

# ── Figure 4: Peak force comparison — all cases at H=0.10 ────────────────
comp = [r for r in results if abs(r["H"]-0.10)<0.001 and r["Fx_peak"] and not np.isnan(r["Fx_peak"])]
if comp:
    labels = [r["label"] for r in comp]
    peaks  = [r["Fx_peak"] for r in comp]
    colors_bar = ['#c55a11' if 'solitary' in r["wtype"] else
                  '#2e75b6' if r["wtype"]=='stokesFirst' else '#70ad47'
                  for r in comp]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(11,5))
    bars = ax.bar(x, peaks, color=colors_bar, edgecolor='black', lw=0.8, width=0.55)
    for bar, v in zip(bars, peaks):
        ax.text(bar.get_x()+bar.get_width()/2, v+2, f'{v:.1f}',
                ha='center', va='bottom', fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Peak Inline Force Fx (N)", fontsize=11)
    ax.set_title("Peak Force Comparison — All Wave Types (H=0.10m, d=0.40m)", fontsize=11)
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#c55a11', label='Solitary'),
                       Patch(facecolor='#2e75b6', label='Stokes 1st'),
                       Patch(facecolor='#70ad47', label='Stokes 2nd')]
    ax.legend(handles=legend_elements, fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    fig.savefig(os.path.join(WORK, "fig4_force_comparison_all.png"), dpi=150, bbox_inches='tight')
    print("Saved: fig4_force_comparison_all.png")

# ── Figure 5: Kt comparison ───────────────────────────────────────────────
kt_cases = [r for r in results if not np.isnan(r["Kt"])]
if kt_cases:
    labels = [r["label"] for r in kt_cases]
    kts    = [r["Kt"] for r in kt_cases]
    colors_kt = ['#c55a11' if 'solitary' in r["wtype"] else
                 '#2e75b6' if r["wtype"]=='stokesFirst' else '#70ad47'
                 for r in kt_cases]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(11,5))
    bars = ax.bar(x, kts, color=colors_kt, edgecolor='black', lw=0.8, width=0.55)
    for bar, v in zip(bars, kts):
        ax.text(bar.get_x()+bar.get_width()/2, v+0.005, f'{v:.3f}',
                ha='center', va='bottom', fontsize=8)
    ax.axhline(1.0, color='red', lw=1, ls='--', label='Kt=1 (no blocking)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Transmission Coefficient Kt", fontsize=11)
    ax.set_title("Wave Transmission Past Cylinder — All Cases", fontsize=11)
    ax.set_ylim(0, 1.3)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    fig.savefig(os.path.join(WORK, "fig5_Kt_comparison.png"), dpi=150, bbox_inches='tight')
    print("Saved: fig5_Kt_comparison.png")


# ── Figure 3b: Peak force vs wave height — solitary vs regular T=3s ───────
reg_H = sorted([r for r in results if r["wtype"]=="stokesFirst" and r["T"]==3.0
                and r["Fx_peak"] and not np.isnan(r["Fx_peak"])], key=lambda r: r["H"])
sol_sorted2 = sorted([r for r in results if r["wtype"]=="solitaryFirst"
                      and r["Fx_peak"] and not np.isnan(r["Fx_peak"])], key=lambda r: r["H"])
if reg_H and sol_sorted2:
    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot([r["H"] for r in sol_sorted2], [r["Fx_peak"] for r in sol_sorted2],
            "o-", color="#c00000", lw=2, ms=9, label="Solitary (solitaryFirst)")
    ax.plot([r["H"] for r in reg_H], [r["Fx_peak"] for r in reg_H],
            "s-", color="#2e75b6", lw=2, ms=9, label="Regular T=3s (stokesFirst)")
    for r in sol_sorted2:
        ax.annotate(f"{r['Fx_peak']:.0f}", (r["H"], r["Fx_peak"]),
                    textcoords="offset points", xytext=(-25,4), fontsize=8, color="#c00000")
    for r in reg_H:
        ax.annotate(f"{r['Fx_peak']:.0f}", (r["H"], r["Fx_peak"]),
                    textcoords="offset points", xytext=(5,4), fontsize=8, color="#2e75b6")
    ax.set_xlabel("Wave Height H (m)", fontsize=11)
    ax.set_ylabel("Peak Inline Force Fx (N)", fontsize=11)
    ax.set_title("Peak Force vs Wave Height — Solitary vs Regular T=3s\n(d=0.40m, cylinder r=0.5m)", fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(WORK, "fig3b_force_vs_H_comparison.png"), dpi=150, bbox_inches="tight")
    print("Saved: fig3b_force_vs_H_comparison.png")


# ── Figure 3b: Peak force vs wave height — solitary vs regular T=3s ───────

reg_H = sorted([r for r in results if r["wtype"]=="stokesFirst" and r["T"]==3.0
                and r["Fx_peak"] and not np.isnan(r["Fx_peak"])], key=lambda r: r["H"])
sol_sorted2 = sorted([r for r in results if r["wtype"]=="solitaryFirst"
                      and r["Fx_peak"] and not np.isnan(r["Fx_peak"])], key=lambda r: r["H"])
if reg_H and sol_sorted2:
    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot([r["H"] for r in sol_sorted2], [r["Fx_peak"] for r in sol_sorted2],
            "o-", color="#c00000", lw=2, ms=9, label="Solitary (solitaryFirst)")
    ax.plot([r["H"] for r in reg_H], [r["Fx_peak"] for r in reg_H],
            "s-", color="#2e75b6", lw=2, ms=9, label="Regular T=3s (stokesFirst)")
    for r in sol_sorted2:
        ax.annotate(f"{r['Fx_peak']:.0f}", (r["H"], r["Fx_peak"]),
                    textcoords="offset points", xytext=(-25,4), fontsize=8, color="#c00000")
    for r in reg_H:
        ax.annotate(f"{r['Fx_peak']:.0f}", (r["H"], r["Fx_peak"]),
                    textcoords="offset points", xytext=(5,4), fontsize=8, color="#2e75b6")
    ax.set_xlabel("Wave Height H (m)", fontsize=11)
    ax.set_ylabel("Peak Inline Force Fx (N)", fontsize=11)
    ax.set_title("Peak Force vs Wave Height — Solitary vs Regular T=3s\n(d=0.40m, cylinder r=0.5m)", fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(WORK, "fig3b_force_vs_H_comparison.png"), dpi=150, bbox_inches="tight")
    print("Saved: fig3b_force_vs_H_comparison.png")


# ── Figure: Peak force vs wave height — Boussinesq vs Chappelear ──────────
chap = sorted([r for r in results if r["wtype"]=="chappelear1962"
               and r["Fx_peak"] and not np.isnan(r["Fx_peak"])], key=lambda r: r["H"])
sol2 = sorted([r for r in results if r["wtype"]=="solitaryFirst"
               and r["Fx_peak"] and not np.isnan(r["Fx_peak"])], key=lambda r: r["H"])
if chap and sol2:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot([r["H"] for r in sol2], [r["Fx_peak"] for r in sol2],
            "o-", color="#c00000", lw=2, ms=9, label="Boussinesq (solitaryFirst)")
    ax.plot([r["H"] for r in chap], [r["Fx_peak"] for r in chap],
            "s--", color="#2e75b6", lw=2, ms=9, label="Chappelear 1962")
    for r in sol2:
        ax.annotate(f"{r['Fx_peak']:.0f}", (r["H"], r["Fx_peak"]),
                    textcoords="offset points", xytext=(-25, 4), fontsize=8, color="#c00000")
    for r in chap:
        ax.annotate(f"{r['Fx_peak']:.0f}", (r["H"], r["Fx_peak"]),
                    textcoords="offset points", xytext=(5, 4), fontsize=8, color="#2e75b6")
    ax.set_xlabel("Wave Height H (m)", fontsize=11)
    ax.set_ylabel("Peak Inline Force Fx (N)", fontsize=11)
    ax.set_title("Solitary Wave Theory Comparison\n(Boussinesq vs Chappelear 1962, d=0.40m)", fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(WORK, "fig_theory_comparison.png"), dpi=150, bbox_inches="tight")
    print("Saved: fig_theory_comparison.png")

# ── save results table ────────────────────────────────────────────────────
with open(os.path.join(WORK, "sweep_results.txt"), 'w') as f:
    f.write(f"{'Case':<22} {'WaveType':<14} {'H':>5} {'T':>5} {'H_inc':>7} {'H_trans':>8} {'Kt':>6} {'Fx_peak':>9} {'Fx_min':>9}\n")
    f.write("-"*85+"\n")
    for r in results:
        T_s = f"{r['T']:.1f}" if r['T'] else "  -"
        Hi_s = f"{r['Hi']:.4f}" if not np.isnan(r['Hi']) else "  nan"
        Ht_s = f"{r['Ht']:.4f}" if not np.isnan(r['Ht']) else "   nan"
        Kt_s = f"{r['Kt']:.3f}" if not np.isnan(r['Kt']) else " nan"
        f.write(f"{r['name']:<22} {r['wtype']:<14} {r['H']:>5.2f} {T_s:>5} "
                f"{Hi_s:>7} {Ht_s:>8} {Kt_s:>6} "
                f"{r['Fx_peak']:>9.2f} {r['Fx_min']:>9.2f}\n")
print("Saved: sweep_results.txt")
