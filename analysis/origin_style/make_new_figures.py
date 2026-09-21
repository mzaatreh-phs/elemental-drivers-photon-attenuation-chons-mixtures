#!/usr/bin/env python3
"""Two new Origin-style figures drawn from figure_data/: (1) carbon-referenced elemental contrast and ensemble sensitivity index versus energy, (2) Geant4-XCOM signed residuals per mixture and mean over the ten mixtures with its standard error.
Also writes their data folders (figure_data/fig_elemental_contrast, fig_validation_residuals) and prints the residual statistics used in the text."""
import sys, os, csv, json
sys.path.insert(0, os.path.dirname(__file__))
from ostyle import *
EL = {'C': (COL[0], 's'), 'H': (COL[1], 'o'), 'N': (COL[2], '^'), 'O': (COL[3], 'v'), 'S': (COL[4], 'D')}
def ser(ax, x, y, color, marker, label=None, ms=5.5, lw=1.0, ls='-', **kw): ax.plot(x, y, color=color, lw=lw, ls=ls, marker=marker, ms=ms * MSF.get(marker, 1.0), mec=color, mfc=color, label=label, **kw)
XL = 'Energy (MeV)'
def writecsv(folder, name, header, rows):
    os.makedirs(FD + folder, exist_ok=True); csv.writer(open(f'{FD}{folder}/{name}.csv', 'w', newline='')).writerows([header] + rows)

# ---------------------------------------------------------------- figure 1: elemental contrast + sensitivity index
h, A = rd('fig_elemental_sensitivity/panelA_MAC_change_per_1wtpct_replacing_C_percent.csv'); h2, I = rd('fig_element_importance/importance_percent.csv')
F1 = 'fig_elemental_contrast'
writecsv(F1, 'panelA_MAC_change_per_percentage_point_replacing_C_percent', ['Energy_MeV', 'H', 'N', 'O', 'S'], [[f'{r[0]:g}'] + ['%.6g' % v for v in r[1:5]] for r in A])
writecsv(F1, 'panelB_sensitivity_index_percent', ['Energy_MeV', 'H', 'N', 'O', 'S'], [[f'{r[0]:g}'] + ['%.4g' % v for v in r[1:5]] for r in I])
open(FD + F1 + '/README.txt', 'w').write("""figs_origin_style/fig_elemental_contrast.png - XCOM mixture rule, base composition = mean of the ten mixtures. Two panels, x log axis (Energy_MeV).
PANEL A (left, symlog y with linear range +-0.3): y = % change of MAC per one percentage point of mass fraction (+1 wt%) of the element replacing carbon -> panelA_MAC_change_per_percentage_point_replacing_C_percent.csv (curves H, N, O, S). Horizontal line at 0.
PANEL B (right, linear y): sensitivity index S_i (%) = |beta_i| / sum_j |beta_j| in the reduced-collinearity ensemble (a normalized absolute standardized sensitivity, NOT a variance share) -> panelB_sensitivity_index_percent.csv (curves H, N, O, S).
""")
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2)); fig.subplots_adjust(wspace=0.27); ax = axs[0]; ax.axhline(0, color='black', lw=0.8, zorder=1)
for k, el in enumerate('HNOS'): c, m = EL[el]; ser(ax, A[:, 0], A[:, 1 + k], c, m, el)
ax.set_yscale('symlog', linthresh=0.3); frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('MAC change per +1 wt% replacing C (%)'); ax.legend(loc='upper right')
ax = axs[1]
for k, el in enumerate('HNOS'): c, m = EL[el]; ser(ax, I[:, 0], I[:, 1 + k], c, m, el)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_ylim(-3, 105); ax.set_xlabel(XL); ax.set_ylabel('Sensitivity index (%)'); ax.legend(loc='upper right')
save(fig, 'fig_elemental_contrast.png')

# ---------------------------------------------------------------- figure 2: validation residuals
h, Dv = rd('fig_validation_deviation/points_deviation_percent.csv'); _, Ev = rd('fig_validation_deviation/errorbars_1sigma_percent.csv'); _, V = rd('fig_validation_deviation/validation_energies_MeV.csv')
E = Dv[:, 0]; dev = Dv[:, 1:]; u = Ev[:, 1:]; n = dev.shape[1]
mean = dev.mean(1); se = np.sqrt((u ** 2).sum(1)) / n; om = dev.mean(); ose = np.sqrt((u ** 2).sum()) / dev.size
F2 = 'fig_validation_residuals'
writecsv(F2, 'panelA_signed_deviation_per_mixture_percent', ['Energy_MeV'] + [f'Mix{m}' for m in range(1, 11)], [[f'{e:g}'] + ['%.6g' % x for x in dev[j]] for j, e in enumerate(E)])
writecsv(F2, 'panelA_errorbars_1sigma_percent', ['Energy_MeV'] + [f'Mix{m}' for m in range(1, 11)], [[f'{e:g}'] + ['%.6g' % x for x in u[j]] for j, e in enumerate(E)])
writecsv(F2, 'panelB_mean_over_mixtures_percent', ['Energy_MeV', 'mean_deviation', 'standard_error_counting'], [[f'{e:g}', '%.6g' % mean[j], '%.6g' % se[j]] for j, e in enumerate(E)])
open(FD + F2 + '/README.txt', 'w').write(f"""figs_origin_style/fig_validation_residuals.png - Geant4 vs XCOM signed MAC deviation. Two panels, x log axis (Energy_MeV).
PANEL A (left): signed deviation (%) of each mixture (10 series) -> panelA_signed_deviation_per_mixture_percent.csv, error bars 1 sigma counting uncertainty -> panelA_errorbars_1sigma_percent.csv; dotted vertical lines = the six validation energies.
PANEL B (right): mean deviation over the ten mixtures at each energy with its standard error sqrt(sum u_i^2)/10 (counting uncertainty only; XCOM uncertainty not included) -> panelB_mean_over_mixtures_percent.csv. Dashed line = overall mean over all 380 points ({om:.3f} %).
""")
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.3)); fig.subplots_adjust(wspace=0.3); ax = axs[0]
for v in V[:, 0]: ax.axvline(v, color='gray', lw=0.8, ls=':', zorder=0)
ax.axhline(0, color='black', lw=0.8, zorder=1)
for i in range(10):
    x = E * (1 + 0.012 * (i - 4.5) / 4.5); ax.errorbar(x, dev[:, i], yerr=u[:, i], fmt='none', ecolor=COL[i], elinewidth=0.8, capsize=1.8, zorder=2); ser(ax, x, dev[:, i], COL[i], MK[i], ms=4, lw=0.0, zorder=3)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('MAC deviation, Geant4 vs XCOM (%)')
ax = axs[1]; ax.axhline(0, color='black', lw=0.8, zorder=1); ax.axhline(om, color=COL[1], lw=1.2, ls='--', zorder=2, label='mean of all points')
ax.errorbar(E, mean, yerr=se, fmt='none', ecolor=COL[2], elinewidth=1.0, capsize=2.5, zorder=2); ser(ax, E, mean, COL[2], 'o', 'mean over ten mixtures', ms=5, lw=0.8, zorder=3)
frame(ax); log10_x(ax, (0.0075, 17)); ax.set_xlabel(XL); ax.set_ylabel('Mean deviation $\\pm$ s.e. (%)'); ax.set_ylim(-0.95, 0.78); ax.legend(loc='lower right')
save(fig, 'fig_validation_residuals.png')
z = mean / se
stats = {'overall_mean': om, 'overall_se': ose, 'mean_min': mean.min(), 'mean_max': mean.max(), 'E_at_min': E[mean.argmin()], 'E_at_max': E[mean.argmax()], 'n_energies_abs_mean_gt_2se': int((abs(z) > 2).sum()), 'n_energies': len(E),
         'n_energies_negative': int((mean < 0).sum()), 'mean_low_E_lt0.05': mean[E < 0.05].mean(), 'mean_E_ge_0.3': mean[E >= 0.3].mean(), 'median_se': float(np.median(se))}
json.dump(stats, open(FD + F2 + '/residual_statistics.json', 'w'), indent=1); print(json.dumps(stats, indent=1))
