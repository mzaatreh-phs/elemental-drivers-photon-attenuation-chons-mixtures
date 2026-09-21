#!/usr/bin/env python3
"""Origin-style versions of the remaining figures (no OriginLab counterpart in figs/): zeff_baseline, process_fractions, validation_deviation, seed_check,
composition_sensitivity, elemental_sensitivity, library_vs_ensemble_correlation, element_importance, exchange_confirmation, broadbeam_*.  Data: figure_data/<fig>/*.csv"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from ostyle import *

XL = 'Energy (MeV)'
EL = {'C': (COL[0], 's'), 'H': (COL[1], 'o'), 'N': (COL[2], '^'), 'O': (COL[3], 'v'), 'S': (COL[4], 'D')}     # element colours / markers
def ser(ax, x, y, color, marker, label=None, ms=5.5, lw=1.0, ls='-', **kw):
    ax.plot(x, y, color=color, lw=lw, ls=ls, marker=marker, ms=ms * MSF.get(marker, 1.0), mec=color, mfc=color, label=label, **kw)
def leg(ax, loc='best', **kw):
    kw.setdefault('borderaxespad', 0.5); return ax.legend(loc=loc, **kw)
def tag(ax, text, loc=(0.03, 0.93)):
    ax.text(loc[0], loc[1], text, transform=ax.transAxes, fontsize=12, va='top', ha='left', bbox=dict(fc='white', ec='none', pad=1.5), zorder=10)
ZERO = dict(color='black', lw=0.8, zorder=1)

# ---------------------------------------------------------------- fig_zeff_baseline
h, Z = rd('fig_zeff_baseline/band_and_reference_Zeff.csv'); E = Z[:, 0]
fig, ax = plt.subplots(figsize=(8.6, 5.7))
ax.fill_between(E, Z[:, 1], Z[:, 2], color=COL[2], alpha=0.25, lw=0, label='CHONS mixtures (range)')
ser(ax, E, Z[:, 3], COL[1], 'o', 'Water (H$_2$O)', ms=6); ser(ax, E, Z[:, 4], COL[3], 's', 'PMMA', ms=6)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('$Z_{eff}$'); ax.set_ylim(3.0, 8.0); leg(ax, 'upper right')
save(fig, 'fig_zeff_baseline.png')

# ---------------------------------------------------------------- fig_process_fractions
h, P = rd('fig_process_fractions/Mix1_process_fractions_percent.csv'); E = P[:, 0]
_, cx = rdtxt('fig_process_fractions/crossover_keV_all_mixtures.csv'); cxk = float([r for r in cx if r[0] == 'Mix1'][0][1])
fig, ax = plt.subplots(figsize=(8.6, 5.7))
for k, (lab, c, m) in enumerate((('Photoelectric', COL[1], 'o'), ('Incoherent (Compton)', COL[2], 's'), ('Coherent', COL[3], '^'), ('Pair production', COL[4], 'D'))): ser(ax, E, P[:, 1 + k], c, m, lab)
ax.axvline(cxk / 1000, color='black', ls='--', lw=1.0, zorder=1); ax.text(cxk / 1000 * 1.08, 93, f'{cxk:.1f} keV', va='center', ha='left', fontsize=12)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_ylim(-4, 104); ax.set_xlabel(XL); ax.set_ylabel('Fraction of total MAC (%)'); leg(ax, 'center', bbox_to_anchor=(0.62, 0.45), borderaxespad=0)
save(fig, 'fig_process_fractions.png')

# ---------------------------------------------------------------- fig_validation_deviation
h, Dv = rd('fig_validation_deviation/points_deviation_percent.csv'); _, Ev = rd('fig_validation_deviation/errorbars_1sigma_percent.csv'); _, V = rd('fig_validation_deviation/validation_energies_MeV.csv')
E = Dv[:, 0]; fig, ax = plt.subplots(figsize=(9.4, 5.7))
for v in V[:, 0]: ax.axvline(v, color='gray', lw=0.8, ls=':', zorder=0)
ax.axhline(0, **ZERO)
for i in range(10):
    x = E * (1 + 0.012 * (i - 4.5) / 4.5)            # tiny horizontal offset so that the ten series do not overlap (data are at the exact energies)
    ax.errorbar(x, Dv[:, 1 + i], yerr=Ev[:, 1 + i], fmt='none', ecolor=COL[i], elinewidth=0.9, capsize=2, zorder=2)
    ser(ax, x, Dv[:, 1 + i], COL[i], MK[i], ms=5, lw=0.0, zorder=3)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('MAC deviation, Geant4 vs XCOM (%)')
ax.legend(handles=mix_handles(), labels=LAB, loc='center left', bbox_to_anchor=(1.02, 0.5))
save(fig, 'fig_validation_deviation.png')

# ---------------------------------------------------------------- fig_seed_check
h, S = rd('fig_seed_check/histogram_and_N01_curve.csv'); w = S[:, 1] - S[:, 0]
fig, ax = plt.subplots(figsize=(6.6, 5.2))
ax.bar(S[:, 2], S[:, 3], width=w, color=COL[2], alpha=0.55, ec='black', lw=0.8, label='seed 1 vs seed 2')
xs = np.linspace(-4.2, 4.2, 300); ax.plot(xs, np.exp(-xs ** 2 / 2) / np.sqrt(2 * np.pi), color=COL[1], lw=2.0, label='N(0,1)')
frame(ax); ax.set_xlim(-4.3, 4.3); ax.set_xlabel('(MAC$_1$ \u2212 MAC$_2$) / \u221a(u$_1^2$ + u$_2^2$)'); ax.set_ylabel('Density'); leg(ax, 'upper right')
save(fig, 'fig_seed_check.png')

# ---------------------------------------------------------------- fig_composition_sensitivity
h, A = rd('fig_composition_sensitivity/panelA_MAC_spread_percent.csv'); h, B = rd('fig_composition_sensitivity/panelB_Pearson_r_with_MAC.csv')
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2)); fig.subplots_adjust(wspace=0.28)
ax = axs[0]; ax.fill_between(A[:, 0], A[:, 2], A[:, 3], color=COL[2], alpha=0.25, lw=0, label='Geant4 68% noise band')
ser(ax, A[:, 0], A[:, 1], COL[2], 'o', 'Geant4'); ax.plot(A[:, 0], A[:, 4], color='black', ls='--', lw=1.4, label='XCOM (noise-free)')
ax.set_yscale('log'); frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('MAC inter-mixture spread (%)'); leg(ax, 'upper right')
ax = axs[1]; ax.axhline(0, **ZERO)
for k, (nm, el) in enumerate((('Oxygen', 'O'), ('Sulfur', 'S'), ('Hydrogen', 'H'))):
    c, m = EL[el]; ser(ax, B[:, 0], B[:, 1 + k], c, m, f'{nm}'); ax.plot(B[:, 0], B[:, 4 + k], color=c, ls='--', lw=1.2)
ax.plot([], [], color='black', ls='--', lw=1.2, label='XCOM (dashed)')
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_ylim(-1.05, 1.05); ax.set_xlabel(XL); ax.set_ylabel('Pearson r with MAC'); leg(ax, 'lower right')
save(fig, 'fig_composition_sensitivity.png')

# ---------------------------------------------------------------- fig_elemental_sensitivity
h, A = rd('fig_elemental_sensitivity/panelA_MAC_change_per_1wtpct_replacing_C_percent.csv'); h, B = rd('fig_elemental_sensitivity/panelB_share_of_MAC_percent_library_mean.csv')
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2)); fig.subplots_adjust(wspace=0.28)
ax = axs[0]; ax.axhline(0, **ZERO)
for k, el in enumerate('HNOS'): c, m = EL[el]; ser(ax, A[:, 0], A[:, 1 + k], c, m, el)
ax.set_yscale('symlog', linthresh=0.3); frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('MAC change per +1 wt% (replacing C) (%)'); leg(ax, 'upper right')
ax = axs[1]
for k, el in enumerate('CHNOS'): c, m = EL[el]; ser(ax, B[:, 0], B[:, 1 + k], c, m, el)
ax.set_yscale('log'); frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('Share of MAC (%)'); leg(ax, 'lower left')
save(fig, 'fig_elemental_sensitivity.png')

# ---------------------------------------------------------------- fig_library_vs_ensemble_correlation
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.9)); fig.subplots_adjust(wspace=0.32)
for ax, (f, nm) in zip(axs, (('panelA_oxygen_r_with_MAC', 'Oxygen'), ('panelB_sulfur_r_with_MAC', 'Sulfur'), ('panelC_hydrogen_r_with_MAC', 'Hydrogen'))):
    h, D = rd(f'fig_library_vs_ensemble_correlation/{f}.csv'); ax.axhline(0, **ZERO)
    ser(ax, D[:, 0], D[:, 1], COL[1], 'o', 'ten-mixture library'); ser(ax, D[:, 0], D[:, 2], COL[2], 's', 'decorrelated ensemble')
    frame(ax); log10_x(ax, (0.0075, 17)); ax.set_ylim(-1.05, 1.05); ax.set_xlabel(XL); ax.set_ylabel(f'r ({nm.lower()} fraction, MAC)')
leg(axs[0], 'lower right')
save(fig, 'fig_library_vs_ensemble_correlation.png')

# ---------------------------------------------------------------- fig_element_importance
h, I = rd('fig_element_importance/importance_percent.csv'); fig, ax = plt.subplots(figsize=(7.6, 5.4))
for k, el in enumerate('HNOS'): c, m = EL[el]; ser(ax, I[:, 0], I[:, 1 + k], c, m, el)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('Importance for MAC variability (%)'); leg(ax, 'upper right')
save(fig, 'fig_element_importance.png')

# ---------------------------------------------------------------- fig_exchange_confirmation
h, rows = rdtxt('fig_exchange_confirmation/panelA_G4_vs_XCOM_MAC_change.csv')
VC = {'dO': ('+5 wt% O',) + EL['O'], 'dN': ('+5 wt% N',) + EL['N'], 'dS': ('+1.5 wt% S',) + EL['S'], 'dH': ('+1.5 wt% H',) + EL['H'], 'dALL': ('all four', 'black', 'o')}
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2)); fig.subplots_adjust(wspace=0.27); ax = axs[0]
ax.plot([-2, 26], [-2, 26], color='black', ls='--', lw=1.0, zorder=1)
for v, (lab, c, m) in VC.items():
    r = np.array([[float(x) for x in q[3:6]] for q in rows if q[0] == v]); ax.errorbar(r[:, 0], r[:, 1], yerr=r[:, 2], fmt=m, color=c, ms=5.5 * MSF.get(m, 1), capsize=2, elinewidth=0.9, label=lab)
frame(ax); ax.set_xlabel('Mixture-rule (XCOM) MAC change (%)'); ax.set_ylabel('Geant4 MAC change (%)'); leg(ax, 'upper left')
ax = axs[1]; h, Zb = rd('fig_exchange_confirmation/panelB_z_vs_energy.csv'); ax.axhline(0, **ZERO)
for y in (-2, 2): ax.axhline(y, color='gray', lw=0.8, ls=':', zorder=1)
for k, (v, (lab, c, m)) in enumerate(VC.items()): ser(ax, Zb[:, 0], Zb[:, 1 + k], c, m, lab, lw=0.0)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_ylim(-2.4, 2.4); ax.set_xlabel(XL); ax.set_ylabel('(Geant4 $-$ XCOM) / uncertainty')
save(fig, 'fig_exchange_confirmation.png')

# ---------------------------------------------------------------- fig_broadbeam_primary_vs_total
fig, axs = plt.subplots(2, 2, figsize=(11.5, 8.4), sharex='col'); fig.subplots_adjust(hspace=0.12, wspace=0.26)
for ax, (f, t) in zip(axs.ravel(), (('panelA_equal10cm_primary_only', 'equal thickness 10 cm, primary only'), ('panelB_equal20gcm2_primary_only', 'equal areal mass 20 g/cm$^2$, primary only'),
                                    ('panelC_equal10cm_with_buildup', 'equal thickness 10 cm, with buildup'), ('panelD_equal20gcm2_with_buildup', 'equal areal mass 20 g/cm$^2$, with buildup'))):
    h, D = rd(f'fig_broadbeam_primary_vs_total/{f}_relative_to_mean_percent.csv'); ax.axhline(0, **ZERO)
    plot_mixes(ax, D[:, 0], D[:, 1:], ms=3.8, lw=0.9); frame(ax); log10_x(ax, (0.0105, 17)); ax.set_ylim(-30, 30); tag(ax, t)
for ax in axs[:, 0]: ax.set_ylabel('Relative transmission (%)')
for ax in axs[1]: ax.set_xlabel(XL)
fig_legend(fig, 'center left', (0.905, 0.5)); save(fig, 'fig_broadbeam_primary_vs_total.png')

# ---------------------------------------------------------------- fig_broadbeam_spread
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2)); fig.subplots_adjust(wspace=0.27)
for ax, (f, t) in zip(axs, (('panelA_equal10cm_spread_percent', 'equal thickness 10 cm'), ('panelB_equal20gcm2_spread_percent', 'equal areal mass 20 g/cm$^2$'))):
    h, D = rd(f'fig_broadbeam_spread/{f}.csv'); ser(ax, D[:, 0], D[:, 1], COL[1], 'o', 'primary only'); ser(ax, D[:, 0], D[:, 2], COL[2], 's', 'with buildup')
    ax.set_yscale('log'); frame(ax); log10_x(ax, (0.0105, 17)); ax.set_xlabel(XL); ax.set_ylabel('Best-to-worst spread (%)'); tag(ax, t, (0.03, 0.06)); leg(ax, 'upper right')
save(fig, 'fig_broadbeam_spread.png')
