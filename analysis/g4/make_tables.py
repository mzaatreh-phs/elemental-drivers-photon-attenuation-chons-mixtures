#!/usr/bin/env python3
"""Build the LaTeX for the spread table (tbl:macspread), correlation table (tbl:correlation) and Zeff table (tbl:zeff) from the analysis outputs
and patch them into ~/CHON/CHOM_submitted version/CHON_shielding(20).tex (Tables 1-3 are handled separately: 1-2 verified, 3 from results/table3_validation.tex)."""
import os, sys, csv, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from g4analysis import *
RES = os.path.dirname(__file__) + '/results/'
TEX = os.path.expanduser('~/CHON/CHOM_submitted version/CHON_shielding(20).tex')
def sg(x): return f'${x:+.3f}$'.replace('-', '-')          # signed coefficient in math mode, as in the original tables

# ------------------------------------------------------------ Table: MAC spread
T4 = json.load(open(RES + 'table4_spreads.json'))
def elab(e): return {0.01: '0.01', 0.03: '0.03', 0.1: '0.10', 1.0: '1.0', 10.0: '10', 15.0: '15'}[e]
def mixid(m, p): return f'{m}' + ('' if p >= 0.9 else r'$^{\dagger}$')
rows = []
for t in T4:
    a, l = t['MAC'], t['LAC']
    rows.append(f"{elab(t['E'])} & {mixid(a['max_mix'], a['p_max_identity'])} & {mixid(a['min_mix'], a['p_min_identity'])} & {a['spread_pct']:.2f} & {a['mc_16']:.2f}--{a['mc_84']:.2f} & {a['xcom_spread_pct']:.2f} & {l['spread_pct']:.1f} \\\\")
tab4 = r'''\begin{table}[pos=H]
\centering
\caption{Selected-energy MAC extrema and relative inter-mixture spread (Eq.~\eqref{eq:spread}) from the 38-energy Geant4 grid, with the 68\% Monte Carlo noise band of the MAC spread (resampling within the counting uncertainty), the noise-free XCOM spread, and the LAC spread. $^{\dagger}$Mixture identity not statistically established (probability below 90\%).}
\label{tbl:macspread}
\scalebox{0.85}{%
\begin{tabular}{ccccccc}
\toprule
 & \multicolumn{2}{c}{MAC extremes} & \multicolumn{3}{c}{MAC spread (\%)} & LAC spread \\
\cmidrule(lr){2-3}\cmidrule(lr){4-6}
Photon energy (MeV) & Max & Min & Geant4 & 68\% noise band & XCOM & (\%) \\
\midrule
''' + '\n'.join(rows) + r'''
\bottomrule
\end{tabular}}
\end{table}'''

# ------------------------------------------------------------ Table: correlations
C = list(csv.DictReader(open(RES + 'table5_correlations.csv')))
EN = [0.01, 0.03, 0.05, 0.08, 0.1, 1.0, 10.0, 15.0]
def cget(pred, resp, key='r_G4'): return [float([r for r in C if r['pred'] == pred and r['resp'] == resp and abs(float(r['E']) - e) < 1e-9][0][key]) for e in EN]
runs, _ = load_runs(); energies = complete_energies(runs); D = build_tables(runs, energies); Eall = list(D['E'])
rho = [RHO[m] for m in MIXES]; ix = [Eall.index(e) for e in EN]
r_mfp = [pearson(rho, list(1 / D['lac'][:, j])) for j in ix]; r_rpe = [pearson(rho, list(100 * (1 - D['T'][:, j]))) for j in ix]
def line(label, vals): return f"{label} & " + ' & '.join(sg(v) for v in vals) + r' \\'
body = '\n'.join([
    line(r'Oxygen fraction vs.\ MAC', cget('Oxygen', 'MAC')), line(r'Hydrogen fraction vs.\ MAC', cget('Hydrogen', 'MAC')),
    line(r'Carbon fraction vs.\ MAC', cget('Carbon', 'MAC')), line(r'Nitrogen (= sulfur = BSA) fraction vs.\ MAC', cget('Nitrogen', 'MAC')), r'\addlinespace',
    line(r'Density vs.\ MAC', cget('Density', 'MAC')), line(r'Density vs.\ LAC', cget('Density', 'LAC')),
    line(r'Density vs.\ MFP, HVL, TVL', r_mfp), line(r'Density vs.\ RPE', r_rpe), r'\addlinespace',
    line(r'Oxygen fraction vs.\ $Z_{\mathrm{eff}}$', cget('Oxygen', 'Zeff')), line(r'Sulfur (= nitrogen = BSA) fraction vs.\ $Z_{\mathrm{eff}}$', cget('Sulfur', 'Zeff')),
    r'\midrule', r'\multicolumn{9}{l}{\emph{Noise-free XCOM values, MAC}} \\',
    line(r'Oxygen fraction vs.\ MAC', cget('Oxygen', 'MAC', 'r_xcom')), line(r'Hydrogen fraction vs.\ MAC', cget('Hydrogen', 'MAC', 'r_xcom')),
    line(r'Carbon fraction vs.\ MAC', cget('Carbon', 'MAC', 'r_xcom')), line(r'Nitrogen (= sulfur = BSA) fraction vs.\ MAC', cget('Nitrogen', 'MAC', 'r_xcom'))])
tab5 = r'''\begin{table}[pos=H]
\centering
\caption{Pearson correlation coefficients ($n=10$ mixtures) between composition or density and the calculated shielding parameters at selected photon energies (Geant4, except the lower block). Nitrogen, sulfur, and BSA mass fractions are perfectly collinear in this library, so one row stands for all three. MFP, HVL, and TVL share one row because they differ only by a constant factor. Above 1~MeV the MAC spread is comparable to the counting noise (Table~\ref{tbl:macspread}); the noise-free XCOM coefficients in the lower block should be used there.}
\label{tbl:correlation}
\footnotesize
\setlength{\tabcolsep}{4pt}
\scalebox{0.85}{%
\begin{tabular}{l cccccccc}
\toprule
 & \multicolumn{8}{c}{Photon energy (MeV)} \\
\cmidrule(lr){2-9}
Variable pair & 0.01 & 0.03 & 0.05 & 0.08 & 0.1 & 1 & 10 & 15 \\
\midrule
''' + body + r'''
\bottomrule
\end{tabular}}
\end{table}'''

# ------------------------------------------------------------ Table: Zeff (XCOM direct method, exact CHONS compositions)
Z = {(r['material'], round(float(r['E_MeV']), 6)): float(r['Zeff_direct']) for r in csv.DictReader(open(os.path.dirname(__file__) + '/../xcom/zeff.csv')) if r['set'] == 'new'}
ZE = [0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.3, 0.6, 1, 2, 5, 10, 15]
zrows = '\n'.join(f"{e:<5g} & " + ' & '.join(f"{Z[(f'Mix{m}', round(e, 6))]:.3f}" for m in MIXES) + r' \\' for e in ZE)
tab6 = r'''\begin{table}[htbp]
\centering
\caption{Effective atomic number of the ten mixtures at selected photon energies (direct method, exact elemental compositions of Table~\ref{tbl:elemental}).}
\label{tbl:zeff}
\small
\setlength{\tabcolsep}{4pt}
\scalebox{0.8}{%
\begin{tabular}{c cccccccccc}
\toprule
 & \multicolumn{10}{c}{$Z_{\mathrm{eff}}$} \\
\cmidrule{2-11}
$E$ (MeV) & Mix~1 & Mix~2 & Mix~3 & Mix~4 & Mix~5 & Mix~6 & Mix~7 & Mix~8 & Mix~9 & Mix~10 \\
\midrule
''' + zrows + r'''
\bottomrule
\end{tabular}}
\end{table}'''
for n, t in (('table4_new.tex', tab4), ('table5_new.tex', tab5), ('table6_zeff_new.tex', tab6)): open(RES + n, 'w').write(t + '\n')

# ------------------------------------------------------------ patch draft (20)
s = open(TEX).read()
def swap(label, new, env_start=r'\begin{table}'):
    global s
    i = s.index(label); a = s.rindex(env_start, 0, i); b = s.index(r'\end{table}', i) + len(r'\end{table}')
    s = s[:a] + new + s[b:]
swap(r'\label{tbl:macspread}', tab4); swap(r'\label{tbl:correlation}', tab5); swap(r'\label{tbl:zeff}', tab6)
open(TEX, 'w').write(s); print('patched', TEX)
