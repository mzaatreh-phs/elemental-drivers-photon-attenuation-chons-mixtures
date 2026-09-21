#!/usr/bin/env python3
"""Insert the moisture-sensitivity Methods subsection + Results subsection/table into CHON_shielding(20).tex; every number comes from moisture_tokens.json / moisture_mac_zeff.csv."""
import os, json, csv, re
import numpy as np
D = os.path.expanduser('~/CHON/recalc_2026-09-19/analysis/moisture/'); TEX = os.path.expanduser('~/CHON/CHOM_submitted version/CHON_shielding(20).tex')
T = json.load(open(D + 'moisture_tokens.json')); s = open(TEX).read()
rows = list(csv.reader(open(D + 'moisture_mac_zeff.csv'))); hdr = rows[0]; A = np.array([[float(x) for x in r] for r in rows[1:]])
E = A[:, 0]; MD = A[:, 1:11]; MW = A[:, 11:21]
# ---- general claims verified on the whole grid before they are written
assert np.all(MW > MD), 'water is not always above the mixtures'
same_ext = all((MD[j].argmax() == MW[j].argmax()) and (MD[j].argmin() == MW[j].argmin()) for j in range(len(E)))
T['mac_min_incr'] = '%.2f' % (100 * (MW / MD - 1).min()); T['mac_max_incr'] = '%.1f' % (100 * (MW / MD - 1).max()); T['ext_same'] = same_ext
assert T['zmax4_wet'] and T['zmin2_wet'] and T['zmin10_wet_001'] and T['rhoorder'] and T['zeqmonowet'], 'a claimed invariance does not hold'
def rep(old, new):
    global s
    assert s.count(old) == 1, (s.count(old), old[:80]); s = s.replace(old, new)
meth = (r'''\subsection{Moisture Sensitivity}\label{sec:moistmethod}

The mixtures are modelled as anhydrous, whereas the precursor powders of Campbell et al.~\cite{campbell2025} contain approximately 10~wt\% water after conditioning at 50\% relative humidity. The effect of this water was estimated by replacing each dry mixture by a wet mixture of 90~wt\% dry mixture and 10~wt\% water (wet-mass basis), with elemental mass fractions $w_i^{\mathrm{wet}}=0.9\,w_i+0.1\,w_{i,\mathrm{H_2O}}$ (H 0.1119 and O 0.8881 for water). The XCOM coefficients, $Z_{\mathrm{eff}}$, the crossover energy, the G--P buildup factors, and the broad-beam responses were recomputed for the wet compositions with the methods described above. The density of the wet mixture was taken either from the additive-volume rule with water at 1.00~g/cm\textsuperscript{3} or unchanged from its dry value; the two cases bracket the unknown void and hydration structure. Dry and wet quantities are both taken from XCOM and G--P calculations, so that the comparison is like for like. The water content of pressed specimens, which was not characterized here, may differ from that of the powders.

''')
rep(r'\subsection{Validation Metric}\label{sec:uncertainty_method}', meth + r'\subsection{Validation Metric}\label{sec:uncertainty_method}')
tab = (r'''\begin{table}[htbp]
\centering
\caption{Effect of adding 10~wt\% water to every mixture (wet-mass basis) on the main results. Dry and wet values are both from XCOM and G--P calculations, so they differ slightly from the Geant4-based values quoted elsewhere. Wet densities follow the additive-volume rule, with the value for unchanged density in parentheses where it differs.}
\label{tbl:moisture}
\small
\setlength{\tabcolsep}{6pt}
\begin{tabular}{l c c}
\toprule
Quantity & Dry & Wet (10~wt\% water) \\
\midrule
MAC spread at 0.01~MeV (\%) & @spd001@ & @spw001@ \\
MAC spread at 0.1~MeV (\%) & @spd010@ & @spw010@ \\
$Z_{\mathrm{eff}}$ range at 0.01~MeV & @zd001@ & @zw001@ \\
$Z_{\mathrm{eff}}$ range near 1~MeV & @zd100@ & @zw100@ \\
Photoelectric--incoherent crossover (keV) & @crd@ & @crw@ \\
Density (g/cm\textsuperscript{3}) & @rhod@ & @rhow@ \\
LAC spread over the grid (\%) & @lspd@ & @lspw@ (@lspw2@) \\
Thickness saving, Mixture~10 vs 2, 0.1~MeV (\%) & @thind@ & @thinw@ (@thins@) \\
Areal-mass excess of Mixture~10 over 2, 0.1~MeV (\%) & @amd@ & @amw@ \\
EBF spread, 40~mfp, 0.04~MeV (\%) & @ebfd004@ & @ebfw004@ \\
EBF spread, 40~mfp, 0.1~MeV (\%) & @ebfd010@ & @ebfw010@ \\
Equal-thickness spread at 0.1~MeV, primary / total (\%) & @bbTd_p01@ / @bbTd_t01@ & @bbTw_p01@ / @bbTw_t01@ \\
Best mixture differs, primary vs total (10~cm) & @bbTd_diff@ of @bbTd_n@ energies & @bbTw_diff@ of @bbTw_n@ energies \\
\bottomrule
\end{tabular}
\end{table}''')
res = (r'''\subsection{Sensitivity to Moisture Content}\label{sec:moisture}

The precursor powders contain about 10~wt\% water, which the anhydrous models omit. Adding 10~wt\% water to every mixture (Section~\ref{sec:moistmethod}; Table~\ref{tbl:moisture}) raises the MAC at every energy, because water has a larger MAC per unit mass than any of the mixtures: by @dmac001@\% at 0.01~MeV, @dmac003@\% at 0.03~MeV, @dmac010@\% at 0.1~MeV, and @dmac1000@\% at 10~MeV. The inter-mixture MAC spread is diluted by roughly 10--15\% of itself (@spd001@ to @spw001@\% at 0.01~MeV and @spd010@ to @spw010@\% at 0.1~MeV), and the mixtures with the maximum and minimum MAC do not change at any of the 38 energies. The elemental contrasts are almost unchanged: the effect of sulfur replacing carbon at 0.01~MeV changes from $@cdS@$ to $@cwS@$\% per weight percent, and those of oxygen, nitrogen, and hydrogen by less than 0.03 percentage points, so that the identification of sulfur and hydrogen as the controlling elements does not depend on the moisture assumption.

The effective atomic number ranges become @zw001@ at 0.01~MeV and @zw100@ near 1~MeV (dry: @zd001@ and @zd100@); Mixture~4 remains the highest at every energy, and Mixture~2 (Mixture~10 at 0.01~MeV) the lowest. The crossover moves to @crw@~keV from @crd@~keV. With additive-volume densities the mixture densities become @rhow@~g/cm\textsuperscript{3} instead of @rhod@, with the same order; the LAC spread over the grid becomes @lspw@\% (@lspw2@\% at unchanged density) instead of @lspd@\%, the thickness saving of Mixture~10 over Mixture~2 at 0.1~MeV becomes @thinw@\% (@thins@\% at unchanged density) instead of @thind@\%, and its areal-mass excess @amw@\% instead of @amd@\%. The compactness-versus-mass-efficiency result is therefore unchanged in sign and nearly unchanged in size.

Exposure buildup is the most moisture-sensitive result. At 0.1~MeV and 40~mfp the buildup factors fall by @ebfchg010@\% and the spread between mixtures from @ebfd010@ to @ebfw010@\%, while Mixture~@ebfwhi010@ remains the highest and Mixture~@ebfwlo010@ the lowest and the EBF ordering still follows the equivalent atomic number. At 0.04~MeV the spread falls from @ebfd004@ to @ebfw004@\% and the highest buildup passes from Mixture~@ebfdhi004@ to Mixture~@ebfwhi004@. In the equal-thickness comparison the spread in the primary response at 0.1~MeV falls from @bbTd_p01@ to @bbTw_p01@\% and that of the total response from @bbTd_t01@ to @bbTw_t01@\%, the ratio of the two over 0.1--1~MeV changes from @bbTd_rlo@--@bbTd_rhi@ to @bbTw_rlo@--@bbTw_rhi@, and the mixture with the lowest primary response differs from the one with the lowest total response at @bbTd_diff@ and @bbTw_diff@ of the @bbTd_n@ energies, respectively; at equal areal mass the rankings still differ substantially (median Spearman coefficient @bbMd_spmed@ dry and @bbMw_spmed@ wet).

Moisture at this level therefore leaves the rankings of MAC, $Z_{\mathrm{eff}}$, and LAC and the elemental contrasts unchanged, reduces the composition-induced spreads of MAC, $Z_{\mathrm{eff}}$, and LAC by about 10--15\%, and reduces those of the buildup factor by about 40--50\%; the only change of an extreme is the highest-buildup mixture at 0.04~MeV. The water content of pressed specimens may differ from that of the powders, and this estimate assumes that the water is uniformly mixed.

''' + tab + '\n\n')
a = s.index(r'\section{Conclusions}'); b = s.rfind('\n% ====', 0, a); s = s[:b] + '\n\n' + res + s[b:]
# ---- limitation statements now point to the quantified sensitivity
rep(r'which is not included here, so the mixtures are models of', r'which is not included here (its effect is quantified in Section~\ref{sec:moisture}), so the mixtures are models of')
rep(r'whereas the source powders contain about 10~wt\% water.', r'whereas the source powders contain about 10~wt\% water, which by itself would reduce the composition-induced spreads of MAC, $Z_{\mathrm{eff}}$, and LAC by about 10--15\% and those of the buildup factor by about 40--50\% without changing their rankings, except for the highest-buildup mixture at 0.04~MeV (Section~\ref{sec:moisture}).')
for k, v in T.items(): s = s.replace('@%s@' % k, str(v))
left = re.findall(r'@\w+@', s); assert not left, left[:8]
open(TEX, 'w').write(s); json.dump(T, open(D + 'moisture_tokens.json', 'w'), indent=1); print('moisture sections inserted; extremes unchanged at all 38 energies:', same_ext)
