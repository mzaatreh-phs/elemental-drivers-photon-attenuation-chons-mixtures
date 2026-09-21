#!/usr/bin/env python3
"""Origin-style versions of the five figures that exist as OriginLab plots in figs/ (mac_lac, mfp_hvl_tvl, rpe, ebf, zeff). Data: figure_data/<fig>/*.csv"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from ostyle import *

XL = 'Energy (MeV)'

# ---------------------------------------------------------------- fig_mac_lac  (Origin: log2 x axis 0.02 ... 16.00 every 2 octaves, linear y, zoom insets)
h, A = rd('fig_mac_lac/panelA_MAC_curves.csv'); _, Au = rd('fig_mac_lac/panelA_MAC_errorbars_1sigma.csv')
_, B = rd('fig_mac_lac/panelB_LAC_curves.csv'); _, Bu = rd('fig_mac_lac/panelB_LAC_errorbars_1sigma.csv')
E = A[:, 0]
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.7)); fig.subplots_adjust(wspace=0.27)
for ax, D, U, yl, xz, xt in ((axs[0], A, Au, 'MAC (cm$^2$/g)', (0.0135, 0.033), [0.02, 0.03]), (axs[1], B, Bu, 'LAC (cm$^{-1}$)', (0.0135, 0.045), [0.02, 0.04])):
    plot_mixes(ax, E, D[:, 1:], yerr=U[:, 1:]); frame(ax)
    log2_x(ax, (2 ** -7.1, 2 ** 4.2), -6, 4, 2); ax.set_xlabel(XL); ax.set_ylabel(yl); pad_ylim(ax, D[:, 1:])
    sel = (E >= xz[0]) & (E <= xz[1]); y = D[sel, 1:]; pad = 0.12 * (y.max() - y.min())
    zoom_inset(ax, xz, (y.min() - pad, y.max() + pad), [0.27, 0.42, 0.36, 0.55], lambda a, D=D, U=U: plot_mixes(a, E, D[:, 1:], ms=7, yerr=U[:, 1:]),
               xfix=lambda a, xz=xz, xt=xt: (a.set_xscale('log', base=2), a.set_xlim(*xz),
                                      a.xaxis.set_major_locator(FixedLocator(xt)), a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}')),
                                      a.xaxis.set_minor_locator(FixedLocator([])), a.tick_params(axis='x', labelsize=12)))
    legend(ax, 'upper right')
save(fig, 'fig_mac_lac.png')

# ---------------------------------------------------------------- fig_mfp_hvl_tvl  (Origin: three tall panels, linear x 0-16, legend at the lower right outside)
fig, axs = plt.subplots(1, 3, figsize=(11.5, 6.6)); fig.subplots_adjust(wspace=0.42)
for ax, f, yl in zip(axs, ('panelA_MFP_cm', 'panelB_HVL_cm', 'panelC_TVL_cm'), ('MFP (cm)', 'HVL (cm)', 'TVL (cm)')):
    h, D = rd(f'fig_mfp_hvl_tvl/{f}.csv'); plot_mixes(ax, D[:, 0], D[:, 1:]); frame(ax)
    ax.set_xlim(-0.6, 16.5); ax.set_xticks(range(0, 17, 2)); ax.set_xlabel(XL); ax.set_ylabel(yl); pad_ylim(ax, D[:, 1:])
fig_legend(fig, 'lower left', (0.905, 0.10)); save(fig, 'fig_mfp_hvl_tvl.png')

# ---------------------------------------------------------------- fig_rpe  (Origin: log10 x, y -10..110, zoom inset 1.1-2.3 MeV)
h, R = rd('fig_rpe/curves_attenuation_1cm_percent.csv'); _, Ru = rd('fig_rpe/errorbars_1sigma_percent.csv'); E = R[:, 0]
fig, ax = plt.subplots(figsize=(8.6, 5.7))
plot_mixes(ax, E, R[:, 1:], yerr=Ru[:, 1:]); frame(ax); log10_x(ax, (0.0075, 17)); ax.set_ylim(-9, 109); ax.set_yticks(range(0, 101, 20))
ax.set_xlabel(XL); ax.set_ylabel('RPE (%)')
xz = (1.1, 2.1); sel = (E >= xz[0]) & (E <= xz[1]); y = R[sel, 1:]; pad = 0.25 * (y.max() - y.min())
yz = (y.min() - pad, y.max() + pad)     # narrow energy window + tight y range so that the differences between mixtures are visible
zoom_inset(ax, xz, yz, [0.33, 0.50, 0.34, 0.44], lambda a: plot_mixes(a, E, R[:, 1:], ms=7, yerr=Ru[:, 1:]),
           xfix=lambda a: (a.set_xticks([1.2, 1.5, 2.0]), a.set_xticklabels(['1.2', '1.5', '2.0']), a.xaxis.set_minor_locator(FixedLocator([]))), corners=('tl', 'br'))
legend(ax, 'upper right'); save(fig, 'fig_rpe.png')

# ---------------------------------------------------------------- fig_ebf  (Origin: 3 x 3 grid, 8 panels + legend in the 9th cell, log10 x 0.01-10, linear y)
fig, axs = plt.subplots(3, 3, figsize=(11.5, 11)); fig.subplots_adjust(hspace=0.46, wspace=0.44)
for k, (ax, d) in enumerate(zip(axs.ravel(), (5, 10, 15, 20, 25, 30, 35, 40))):
    h, D = rd(f'fig_ebf/panel{k + 1}_{d}mfp_EBF.csv'); plot_mixes(ax, D[:, 0], D[:, 1:], ms=4.8); frame(ax)
    log10_x(ax, (0.01, 17)); ax.set_xlabel(XL); ax.set_ylabel(f'EBF (at {d} mfp)'); pad_ylim(ax, D[:, 1:], lo_frac=0.05)
    ax.ticklabel_format(axis='y', style='plain')
axs[2, 2].axis('off'); legend(axs[2, 2], 'center', bbox_to_anchor=(0.5, 0.5), borderaxespad=0)
save(fig, 'fig_ebf.png')

# ---------------------------------------------------------------- fig_zeff  (Origin: log2 x, every octave 0.02 ... 16.00, big markers, legend upper right)
h, Z = rd('fig_zeff/curves_Zeff.csv'); E = Z[:, 0]
fig, ax = plt.subplots(figsize=(10, 6.8))
plot_mixes(ax, E, Z[:, 1:], ms=9.5, lw=1.0); frame(ax); log2_x(ax, (2 ** -6.9, 2 ** 4.2), -6, 4, 1)
ax.set_xlabel('Photon Energy (MeV)'); ax.set_ylabel('$Z_{eff}$'); ax.set_ylim(3.35, 7.7); ax.set_yticks(np.arange(3.5, 7.6, 0.5))
legend(ax, 'upper right'); save(fig, 'fig_zeff.png')
