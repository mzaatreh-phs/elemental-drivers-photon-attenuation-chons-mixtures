#!/usr/bin/env python3
"""Figures for the recalculated CHON study -> analysis/g4/figs/.  XCOM-only figures are always made; Geant4 figures use the
energies completed so far (all 10 mixtures x 2 seeds)."""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from g4analysis import *

FIG = f'{OUT}/figs/'; os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.25, 'figure.dpi': 120, 'savefig.dpi': 300,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
COL = plt.get_cmap('tab10').colors; MK = ['o', 's', '^', 'v', 'D', '<', '>', 'p', 'h', '*']
def save(fig, name): fig.savefig(FIG + name, bbox_inches='tight'); plt.close(fig); print('wrote', name)
def lines(ax, x, Y, ms=3, lw=0.9):
    for i, m in enumerate(MIXES): ax.plot(x, Y[i], marker=MK[i], ms=ms, lw=lw, color=COL[i], label=f'Mix {m}')

zeffd, zref = load_zeff(); xc = load_xcom_mac()
runs, _ = load_runs(); energies = complete_energies(runs)

# ---------------- XCOM-only figures ----------------
Eg = np.array(GRID)
Zm = np.array([[zeffd[(m, round(e, 6))] for e in Eg] for m in MIXES])
fig, ax = plt.subplots(figsize=(5.2, 3.6)); lines(ax, Eg, Zm); ax.set_xscale('log')
ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel(r'$Z_\mathrm{eff}$'); ax.legend(ncol=2, fontsize=7); save(fig, 'fig_zeff.png')

lo, hi = Zm.min(0), Zm.max(0)
fig, ax = plt.subplots(figsize=(5.2, 3.6)); ax.fill_between(Eg, lo, hi, alpha=0.3, color='tab:blue', label='CHONS mixtures (range)')
for nm, c in (('water', 'tab:red'), ('PMMA', 'tab:green')):
    ax.plot(Eg, [zref[(nm, round(e, 6))] for e in Eg], marker='o', ms=3, color=c, label={'water': r'Water (H$_2$O)', 'PMMA': 'PMMA'}[nm])
ax.set_xscale('log'); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel(r'$Z_\mathrm{eff}$'); ax.legend(); save(fig, 'fig_zeff_baseline.png')

P = {}
for r in csv.DictReader(open(f'{XC}/process_fractions.csv')):
    if r['set'] == 'new' and r['material'] == 'Mix1': P[round(float(r['E_MeV']), 6)] = r
cx = [float(r['crossover_keV']) for r in csv.DictReader(open(f'{XC}/crossover.csv')) if r['set'] == 'new' and r['material'] == 'Mix1'][0]
fig, ax = plt.subplots(figsize=(5.2, 3.6))
for k, lab in (('F_photoelectric', 'Photoelectric'), ('F_incoherent', 'Incoherent'), ('F_coherent', 'Coherent'), ('F_pair', 'Pair production')):
    ax.plot(Eg, [100 * float(P[round(e, 6)][k]) for e in Eg], marker='o', ms=3, label=lab)
ax.axvline(cx / 1000, ls='--', color='k', lw=0.8); ax.text(cx / 1000 * 1.05, 50, f'{cx:.1f} keV', rotation=90, va='center', fontsize=8)
ax.set_xscale('log'); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel('Fraction of total MAC (%)'); ax.legend(fontsize=7); save(fig, 'fig_process_fractions.png')

EB = {}
for r in csv.DictReader(open(f'{XC}/exposure_buildup.csv')):
    if r['set'] == 'new': EB[(int(r['material'][3:]), round(float(r['E_MeV']), 6))] = r
E37 = [e for e in GRID if e >= 0.015]
fig, axs = plt.subplots(2, 4, figsize=(11, 5.2), sharex=True)
for ax, d in zip(axs.ravel(), (5, 10, 15, 20, 25, 30, 35, 40)):
    for i, m in enumerate(MIXES): ax.plot(E37, [float(EB[(m, round(e, 6))][f'B_{d}mfp']) for e in E37], color=COL[i], lw=0.9, label=f'Mix {m}')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_title(f'{d} mfp', fontsize=9)
for ax in axs[1]: ax.set_xlabel('Energy (MeV)')
for ax in axs[:, 0]: ax.set_ylabel('EBF')
axs[0, 0].legend(fontsize=6, ncol=2); fig.tight_layout(); save(fig, 'fig_ebf.png')

# ---------------- Geant4 figures ----------------
if energies:
    D = build_tables(runs, energies); E = D['E']
    Xc = np.array([[xc[(m, round(e, 6))] for e in E] for m in MIXES])
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    lines(axs[0], E, D['mac']); lines(axs[1], E, D['lac'])
    for ax, yl in zip(axs, ('MAC (cm$^2$/g)', 'LAC (cm$^{-1}$)')): ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel(yl)
    axs[0].legend(ncol=2, fontsize=7); save(fig, 'fig_mac_lac.png')

    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, (f, yl) in zip(axs, ((lambda l: 1 / l, 'MFP (cm)'), (lambda l: np.log(2) / l, 'HVL (cm)'), (lambda l: np.log(10) / l, 'TVL (cm)'))):
        lines(ax, E, f(D['lac'])); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel(yl)
    axs[0].legend(ncol=2, fontsize=6); fig.tight_layout(); save(fig, 'fig_mfp_hvl_tvl.png')

    fig, ax = plt.subplots(figsize=(5.2, 3.6)); lines(ax, E, 100 * (1 - D['T'])); ax.set_xscale('log')
    ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel('Primary-beam attenuation, 1 cm (%)'); ax.legend(ncol=2, fontsize=7); save(fig, 'fig_rpe.png')

    dp = 100 * (D['mac'] / Xc - 1); sd = 100 * D['u'] / Xc
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for i, m in enumerate(MIXES): ax.errorbar(E * (1 + 0.012 * (i - 4.5) / 4.5), dp[:, :][i], yerr=sd[i], fmt=MK[i], ms=3, lw=0.6, color=COL[i], label=f'Mix {m}')
    ax.axhline(0, color='k', lw=0.6); ax.set_xscale('log'); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel('Geant4 vs XCOM, MAC deviation (%)')
    for v in VALID: ax.axvline(v, color='gray', lw=0.4, ls=':')
    ax.legend(ncol=5, fontsize=6, loc='lower right'); save(fig, 'fig_validation_deviation.png')

    z = (D['mac1'] - D['mac2']) / np.sqrt(D['u1'] ** 2 + D['u2'] ** 2)
    fig, ax = plt.subplots(figsize=(4.4, 3.2)); ax.hist(z.ravel(), bins=np.linspace(-4, 4, 33), density=True, alpha=0.6, label='seed 1 vs seed 2')
    xs = np.linspace(-4, 4, 200); ax.plot(xs, np.exp(-xs ** 2 / 2) / np.sqrt(2 * np.pi), 'k', lw=1, label='N(0,1)')
    ax.set_xlabel(r'$(\mu_{m,1}-\mu_{m,2})/\sqrt{u_1^2+u_2^2}$'); ax.set_ylabel('Density'); ax.legend(); save(fig, 'fig_seed_check.png')

    # composition sensitivity vs energy: MAC spread (with noise band) and O/S/H correlation with MAC (G4 and noise-free XCOM)
    sp_m, sp_lo, sp_hi, sp_x = [], [], [], []
    for j in range(len(E)):
        s = np.array([spread_pct(d) for d in mc_draw(D['mac'][:, j], D['u'][:, j], n=1500)]); l, m_, h = interval(s)
        sp_m.append(spread_pct(D['mac'][:, j])); sp_lo.append(l); sp_hi.append(h); sp_x.append(spread_pct(Xc[:, j]))
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.6))
    axs[0].fill_between(E, sp_lo, sp_hi, alpha=0.3); axs[0].plot(E, sp_m, 'o-', ms=3, label='Geant4 (68% noise band)'); axs[0].plot(E, sp_x, 'k--', lw=1, label='XCOM (noise-free)')
    axs[0].set_xscale('log'); axs[0].set_yscale('log'); axs[0].set_xlabel('Photon energy (MeV)'); axs[0].set_ylabel('MAC inter-mixture spread (%)'); axs[0].legend(fontsize=7)
    PR = {'Oxygen': [COMP[m]['O'] for m in MIXES], 'Sulfur': [COMP[m]['S'] for m in MIXES], 'Hydrogen': [COMP[m]['H'] for m in MIXES]}
    for k, (nm, pv) in enumerate(PR.items()):
        rg = [pearson(pv, D['mac'][:, j]) for j in range(len(E))]; rx = [pearson(pv, Xc[:, j]) for j in range(len(E))]
        axs[1].plot(E, rg, 'o-', ms=3, color=COL[k], label=f'{nm} (Geant4)'); axs[1].plot(E, rx, '--', color=COL[k], lw=1)
    axs[1].axhline(0, color='k', lw=0.6); axs[1].set_xscale('log'); axs[1].set_ylim(-1.05, 1.05); axs[1].set_xlabel('Photon energy (MeV)'); axs[1].set_ylabel('Pearson r with MAC'); axs[1].legend(fontsize=7)
    fig.tight_layout(); save(fig, 'fig_composition_sensitivity.png')
print('energies used for Geant4 figures:', len(energies))
