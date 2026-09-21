#!/usr/bin/env python3
import os, csv, json, shutil
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = os.path.expanduser('~/CHON/recalc_2026-09-19'); S = R + '/analysis/sensitivity/'; B = R + '/analysis/broadbeam/'
DEST = os.path.expanduser('~/CHON/CHOM_submitted version/figs_new/'); DATA = os.path.expanduser('~/CHON/CHOM_submitted version/figure_data/novelty_items_1_2/')
# data export moved to analysis/g4/export_figure_data.py (one folder per figure)
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.25, 'savefig.dpi': 300, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
COL = plt.get_cmap('tab10').colors
def rd(fn): 
    r = list(csv.reader(open(fn))); return r[0], np.array([[float(x) if x != '' else np.nan for x in row] for row in r[1:]])
# ---------------- item 1
h, a = rd(S + 'exchange_sensitivity_and_contribution.csv'); E = a[:, 0]
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.6))
for k, (el, c) in enumerate(zip(('H', 'N', 'O', 'S'), COL)): axs[0].plot(E, a[:, 1 + k], 'o-', ms=3, color=c, label=el)
axs[0].axhline(0, color='k', lw=0.6); axs[0].set_xscale('log'); axs[0].set_yscale('symlog', linthresh=0.3); axs[0].set_xlabel('Photon energy (MeV)'); axs[0].set_ylabel('MAC change per +1 wt% (replacing C) (%)'); axs[0].legend()
for k, (el, c) in enumerate(zip(('C', 'H', 'N', 'O', 'S'), COL)): axs[1].plot(E, a[:, 5 + k], 'o-', ms=3, color=c, label=el)
axs[1].set_xscale('log'); axs[1].set_yscale('log'); axs[1].set_xlabel('Photon energy (MeV)'); axs[1].set_ylabel('Share of MAC (%) (library mean)'); axs[1].legend()
fig.tight_layout(); fig.savefig(DEST + 'fig_elemental_sensitivity.png', bbox_inches='tight'); plt.close(fig)
h, c = rd(S + 'ensemble_vs_library_correlation_and_importance.csv'); E = c[:, 0]
fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
for ax, (k, el) in zip(axs, ((2, 'O'), (3, 'S'), (0, 'H'))):
    ax.plot(E, c[:, 1 + 4 * k], 'o-', ms=3, color='tab:red', label='ten-mixture library'); ax.plot(E, c[:, 2 + 4 * k], 's-', ms=3, color='tab:blue', label='decorrelated ensemble')
    ax.axhline(0, color='k', lw=0.6); ax.set_xscale('log'); ax.set_ylim(-1.05, 1.05); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel(f'r({el} fraction, MAC)'); ax.set_title(el, fontsize=9)
axs[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(DEST + 'fig_library_vs_ensemble_correlation.png', bbox_inches='tight'); plt.close(fig)
fig, ax = plt.subplots(figsize=(5.4, 3.6))
for k, (el, col) in enumerate(zip('HNOS', COL)): ax.plot(E, c[:, 4 + 4 * k], 'o-', ms=3, color=col, label=el)
ax.set_xscale('log'); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel('Importance for MAC variability (%)'); ax.legend(ncol=4); fig.savefig(DEST + 'fig_element_importance.png', bbox_inches='tight'); plt.close(fig)
# ---------------- item 2
def wide(fn): 
    hh, aa = rd(fn); return aa[:, 0], aa[:, 1:]
fig, axs = plt.subplots(2, 2, figsize=(10, 6.2), sharex='col')
for j, (tag, ttl) in enumerate((('t10cm', 'equal thickness 10 cm'), ('m20gcm2', 'equal areal mass 20 g/cm$^2$'))):
    Ep, Tp = wide(B + f'wide_T_primary_{tag}.csv'); Eb, Tb = wide(B + f'wide_T_broad_{tag}.csv')
    for row, (Ee, T, lab) in enumerate(((Ep, Tp, 'primary only'), (Eb, Tb, 'with buildup'))):
        ax = axs[row, j]; rel = 100 * (T / np.nanmean(T, axis=1, keepdims=True) - 1)
        for m in range(10): ax.plot(Ee, rel[:, m], color=COL[m], lw=0.9, marker='o', ms=2, label=f'Mix {m + 1}')
        ax.axhline(0, color='k', lw=0.6); ax.set_xscale('log'); ax.set_title(f'{ttl}: {lab}', fontsize=9); ax.set_ylabel('Transmission relative to mixture mean (%)')
        ax.set_ylim(-30, 30)
axs[1, 0].set_xlabel('Photon energy (MeV)'); axs[1, 1].set_xlabel('Photon energy (MeV)'); axs[0, 0].legend(ncol=2, fontsize=6)
fig.tight_layout(); fig.savefig(DEST + 'fig_broadbeam_primary_vs_total.png', bbox_inches='tight'); plt.close(fig)
r = list(csv.DictReader(open(B + 'broadbeam_summary.csv')))
fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.5))
for ax, (kind, v, ttl) in zip(axs, (('thickness_cm', '10', 'equal thickness 10 cm'), ('areal_mass_g_cm2', '20', 'equal areal mass 20 g/cm$^2$'))):
    s = [x for x in r if x['condition'] == kind and x['value'] == v]; e = [float(x['E_MeV']) for x in s]
    ax.plot(e, [float(x['spread_primary_pct']) for x in s], 'o-', ms=3, label='primary only'); ax.plot(e, [float(x['spread_broad_pct']) for x in s], 's-', ms=3, label='with buildup')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Photon energy (MeV)'); ax.set_ylabel('Best-to-worst mixture spread in transmission (%)'); ax.set_title(ttl, fontsize=9); ax.legend()
fig.tight_layout(); fig.savefig(DEST + 'fig_broadbeam_spread.png', bbox_inches='tight'); plt.close(fig)
print('figures ->', DEST)
