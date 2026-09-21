#!/usr/bin/env python3
"""Insert the scenario-specification table (Methods), the dominant-element table and the sigma_O/sigma_S criterion (Results) into CHON_shielding(20).tex."""
import os, csv, json, re
import numpy as np
A = os.path.expanduser('~/CHON/recalc_2026-09-19/analysis/sensitivity/'); TEX = os.path.expanduser('~/CHON/CHOM_submitted version/CHON_shielding(20).tex')
s = open(TEX).read()
sc = {r['scenario'][0]: r for r in csv.DictReader(open(A + 'scenario_correlations.csv'))}
rows = list(csv.DictReader(open(A + 'element_mac.csv'))); E = np.array([float(r['Energy_MeV']) for r in rows]); EL = ['C', 'H', 'N', 'O', 'S']
MU = np.array([[float(r[f'mu_{e}_cm2_per_g']) for r in rows] for e in EL]); j = int(np.argmin(abs(E - 0.01))); crit = abs((MU[4, j] - MU[0, j]) / (MU[3, j] - MU[0, j]))
inp = json.load(open(os.path.expanduser('~/CHON/recalc_2026-09-19/inputs/inputs_exact.json')))['mixtures']
W = np.array([[inp[str(m)]['composition'].get(e, 0.0) for e in EL] for m in range(1, 11)]); lib_ratio = W[:, 3].std(ddof=1) / W[:, 4].std(ddof=1)
ratio = {k: float(v['sd_O']) / float(v['sd_S']) for k, v in sc.items()}
def rep(old, new):
    global s
    assert s.count(old) == 1, (s.count(old), old[:80]); s = s.replace(old, new)
# ---- Methods: shorten the prose enumeration, refer to the table
a = s.index(r'Because such rankings depend on the assumed compositions, they were repeated for'); b = s.index(r'In scenarios A--F the dominant element was identified')
s = s[:a] + r'''Because such rankings depend on the assumed compositions, they were repeated for seven scenarios (A--G; Table~\ref{tbl:scenarios}): six sampled with $2\times10^{4}$ compositions each and one analytic (G: independent inputs of equal variance, for which the variance shares are proportional to $\Delta_i^{2}$). The scenarios differ in the spread of each element, in the condition placed on carbon, and hence in the residual correlations between the elements, which are listed in the table; the standard deviations and correlations of every scenario are provided with the data repository. ''' + s[b:]
tab1 = (r'''\begin{table}[htbp]
\centering
\caption{Scenarios used to test the dependence of the element ranking on the assumed compositions. The last two columns are computed from the sampled compositions: the ratio of the standard deviations of the oxygen and sulfur mass fractions, and the largest absolute correlation between two of H, N, O, and S. Scenario A is the ensemble described in the text regenerated with $2\times10^{4}$ compositions.}
\label{tbl:scenarios}
\footnotesize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{l p{5.6cm} p{3.4cm} c c}
\toprule
Scenario & Sampling of H, N, O, and S & Condition on carbon (realized range) & $\sigma_{\mathrm{O}}/\sigma_{\mathrm{S}}$ & max $|r|$ \\
\midrule
A & uniform between the library minimum and maximum & within library range (@Alo@--@Ahi@) & @rA@ & @mA@ \\
B & uniform, half-ranges halved about the library mean (clipped at 0) & $>0.2$ (@Blo@--@Bhi@) & @rB@ & @mB@ \\
C & uniform, half-ranges doubled about the library mean (clipped at 0) & $>0.2$ (@Clo@--@Chi@) & @rC@ & @mC@ \\
D & as A without the carbon condition, with S $=0.114$ N (BSA ratio) & $>0.2$ (@Dlo@--@Dhi@) & @rD@ & @mD@ \\
E & Dirichlet around the library mean, concentration 300 & none (closure; @Elo@--@Ehi@) & @rE@ & @mE@ \\
F & uniform: H 0.04--0.14, N 0--0.25, O 0--0.60, S 0--0.05 & 0.2--0.85 (@Flo@--@Fhi@) & @rF@ & @mF@ \\
G & analytic: independent inputs of equal variance & --- & 1 & 0 \\
\bottomrule
\end{tabular}
\end{table}''')
for k in 'ABCDEF':
    r = sc[k]; tab1 = tab1.replace('@%slo@' % k, '%.2f' % float(r['carbon_min'])).replace('@%shi@' % k, '%.2f' % float(r['carbon_max'])).replace('@r%s@' % k, '%.1f' % ratio[k]).replace('@m%s@' % k, '%.2f' % float(r['max_abs_offdiag_r']))
a = s.index(r'In scenarios A--F the dominant element was identified'); b = s.index('\n\n', a); s = s[:b] + '\n\n' + tab1 + s[b:]
# ---- Results: criterion paragraph + dominant-element table
R = list(csv.DictReader(open(A + 'range_test_importance.csv'))); TE = (0.01, 0.03, 0.05, 0.1, 1, 3, 10, 15); names = {'A': 'A', 'B': 'B', 'C': 'C', 'D': 'D', 'E': 'E', 'F': 'F', 'G': 'G'}
def dom(sname, e):
    r = [x for x in R if x['scenario'] == sname and abs(float(x['Energy_MeV']) - e) < 1e-9][0]; sh = {el: float(r['share_' + el]) for el in 'HNOS'}; return max(sh, key=sh.get)
scn = sorted({r['scenario'] for r in R}); table_rows = '\n'.join('%s & ' % sn[0] + ' & '.join(dom(sn, e) for e in TE) + r' \\' for sn in scn)
para = (r'''The dependence on the scenario can be stated quantitatively. For nearly independent inputs the variance of the MAC at a given energy is proportional to $\Delta_i^{2}\sigma_i^{2}$, so at 0.01~MeV, where the contrast of sulfur per unit mass is @crit@ times that of oxygen, sulfur dominates whenever the spread of the oxygen mass fraction is less than @crit@ times that of sulfur. The ratio $\sigma_{\mathrm{O}}/\sigma_{\mathrm{S}}$ is @rA@ in scenario A, @rE@ in E, and @rF@ in F, where sulfur dominates, and @rB@ in B and @rC@ in C, where oxygen dominates (Table~\ref{tbl:scenarios}); the ten recipes themselves have a ratio of @lib@, close to this boundary. Sulfur dominance of the MAC variability is therefore a statement about the chosen compositional spreads, whereas its larger contrast per unit mass is not (Table~\ref{tbl:dominance}).

\begin{table}[htbp]
\centering
\caption{Element with the largest variance share of the MAC (H, N, O, or S) in each scenario of Table~\ref{tbl:scenarios} at selected photon energies. Scenarios A--F use the variance share $\Delta_i\,\mathrm{cov}(w_i,\mu_m^{\mathrm{mix}})/\mathrm{var}(\mu_m^{\mathrm{mix}})$; G is analytic.}
\label{tbl:dominance}
\small
\setlength{\tabcolsep}{7pt}
\begin{tabular}{c cccccccc}
\toprule
 & \multicolumn{8}{c}{Photon energy (MeV)} \\
\cmidrule(lr){2-9}
Scenario & 0.01 & 0.03 & 0.05 & 0.1 & 1 & 3 & 10 & 15 \\
\midrule
''' + table_rows + r'''
\bottomrule
\end{tabular}
\end{table}''')
para = para.replace('@crit@', '%.1f' % crit).replace('@lib@', '%.1f' % lib_ratio)
for k in 'ABCEF': para = para.replace('@r%s@' % k, '%.1f' % ratio[k])
end = r'and that oxygen is not a general control variable.'
a = s.index(end); b = s.index('\n\n', a); s = s[:b] + '\n\n' + para + s[b:]
left = re.findall(r'@\w+@', s); assert not left, left
open(TEX, 'w').write(s); print('scenario tables inserted; crit=%.1f lib=%.1f' % (crit, lib_ratio))
