#!/usr/bin/env python3
"""Numbers for Results 3.2 (sensitivity, process), 3.6 (Zeff), 3.7 (EBF, broad-beam)."""
import os, sys, csv, json
import numpy as np
A = os.path.expanduser('~/CHON/recalc_2026-09-19/analysis/')
def rows(p): return list(csv.DictReader(open(A + p)))
print('== LAC largest not Mix10 at energies:')
sys.path.insert(0, A + 'g4')
from g4analysis import *
runs, _ = load_runs(); en = complete_energies(runs); Dd = build_tables(runs, en); E = np.array(Dd['E'])
print('  ', [(float(E[j]), MIXES[Dd['lac'][:, j].argmax()]) for j in range(len(E)) if Dd['lac'][:, j].argmax() != 9])
print('== elemental sensitivity: % MAC change per +1 wt% replacing C')
S = rows('sensitivity/exchange_sensitivity_and_contribution.csv'); k = list(S[0].keys())
for r in S:
    if float(r['Energy_MeV']) in (0.01, 0.02, 0.03, 0.05, 0.1, 1, 3, 10, 15): print('  E=%-5s' % r['Energy_MeV'], {kk.split('_')[4]: round(float(v), 3) for kk, v in r.items() if kk.startswith('pct_MAC')})
print('== ensemble vs library correlation and importance (H,N,O,S)')
C = rows('sensitivity/ensemble_vs_library_correlation_and_importance.csv')
for r in C:
    if float(r['Energy_MeV']) in (0.01, 0.03, 0.05, 0.1, 1, 10, 15): print('  E=%-5s' % r['Energy_MeV'], 'r_lib O %+.2f S %+.2f H %+.2f | r_ens O %+.2f S %+.2f H %+.2f | imp%% H %.0f N %.0f O %.0f S %.0f' % tuple(float(r[x]) for x in ('r_library_O', 'r_library_S', 'r_library_H', 'r_ensemble_O', 'r_ensemble_S', 'r_ensemble_H', 'importance_pct_H', 'importance_pct_N', 'importance_pct_O', 'importance_pct_S')))
print('== range test: scenarios and where the dominant element (max share) is H / S / O')
R = rows('sensitivity/range_test_importance.csv'); sc = sorted({r['scenario'] for r in R}); print('  scenarios:', sc)
for e in (0.01, 0.03, 0.05, 0.1, 1, 3, 10, 15):
    dom = {}
    for s_ in sc:
        r = [x for x in R if x['scenario'] == s_ and abs(float(x['Energy_MeV']) - e) < 1e-9]
        if r: r = r[0]; sh = {el: float(r['share_' + el]) for el in 'HNOS'}; dom[s_[:1]] = max(sh, key=sh.get)
    print('   E=%-5g dominant by scenario:' % e, dom)
print('== exchange confirmation (Geant4 vs mixture rule): examples and z stats')
X = [r for r in rows('sensitivity/exchange_confirmation.csv') if r['variant'] in ('dO', 'dN', 'dS', 'dH', 'dALL')]
z = np.array([float(r['z']) for r in X]); print('  n=%d mean z %.2f std z %.2f chi2/N %.2f max|z| %.2f' % (len(z), z.mean(), z.std(ddof=1), (z ** 2).mean(), abs(z).max()))
for r in X:
    if float(r['E_MeV']) in (0.01, 0.03, 0.1) and r['variant'] in ('dS', 'dO', 'dALL', 'dH', 'dN'): print('   E=%-5s %-5s G4 %+.2f +-%.2f  XCOM %+.2f  z %+.2f' % (r['E_MeV'], r['variant'], float(r['dMAC_Geant4_pct']), float(r['u_pct']), float(r['dMAC_XCOM_mixture_rule_pct']), float(r['z'])))
print('== process fractions Mix1')
P = [r for r in rows('xcom/process_fractions.csv') if r['set'] == 'new' and r['material'] == 'Mix1']
for r in P:
    if float(r['E_MeV']) in (0.01, 0.02, 0.03, 0.1, 1, 10, 15): print('   E=%-5s PE %.1f%% INC %.1f%% COH %.1f%% PAIR %.1f%%' % (r['E_MeV'], 100 * float(r['F_photoelectric']), 100 * float(r['F_incoherent']), 100 * float(r['F_coherent']), 100 * float(r['F_pair'])))
coh = [(float(r['E_MeV']), 100 * float(r['F_coherent'])) for r in P]; print('   coherent peak (Mix1):', max(coh, key=lambda t: t[1]))
CR = {r['material']: float(r['crossover_keV']) for r in rows('xcom/crossover.csv') if r['set'] == 'new'}; print('   crossover keV range %.2f-%.2f ; Mix1 %.2f ; %s' % (min(CR.values()), max(CR.values()), CR['Mix1'], {k: round(v, 2) for k, v in CR.items()}))
cohmax = max(100 * float(r['F_coherent']) for r in rows('xcom/process_fractions.csv') if r['set'] == 'new'); print('   max coherent fraction over all mixtures %.1f%%' % cohmax)
print('== Zeff: per-energy range, extremes, water/PMMA vs band')
Z = {}
for r in rows('xcom/zeff.csv'):
    if r['set'] in ('new', 'ref'): Z[(r['material'], float(r['E_MeV']))] = float(r['Zeff_direct'])
Es = sorted({k[1] for k in Z})
for e in Es:
    v = np.array([Z[('Mix%d' % m, e)] for m in range(1, 11)]); w, p = Z[('water', e)], Z[('PMMA', e)]
    pos = lambda x: 'above' if x > v.max() else ('below' if x < v.min() else 'within')
    if e in (0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.1, 0.3, 1, 2, 5, 10, 15): print('   E=%-5g band %.3f (Mix%d)-%.3f (Mix%d) spread %.3f | water %.3f %s | PMMA %.3f %s' % (e, v.min(), v.argmin() + 1, v.max(), v.argmax() + 1, v.max() - v.min(), w, pos(w), p, pos(p)))
ch = []
last = None
for e in Es:
    v = np.array([Z[('Mix%d' % m, e)] for m in range(1, 11)]); o = (v.argmax() + 1, v.argmin() + 1)
    if o != last: ch.append((e, o)); last = o
print('   (max,min) mixture changes:', ch)
print('   Zeff Mix4 vs Mix10:', [(e, round(Z[('Mix4', e)], 3), round(Z[('Mix10', e)], 3)) for e in (0.03, 0.04, 0.05, 0.08, 1, 5, 10, 15)])
print('   water/PMMA classification changes:', [(e, ('above' if Z[('water', e)] > max(Z[('Mix%d' % m, e)] for m in range(1, 11)) else 'below' if Z[('water', e)] < min(Z[('Mix%d' % m, e)] for m in range(1, 11)) else 'within')) for e in Es][::3])
inpj = json.load(open(A + '../inputs/inputs_exact.json'))['mixtures']; print('   mean hydrogen wt%% of mixtures %.1f (range %.1f-%.1f); water H 11.19' % (100 * np.mean([inpj[str(m)]['composition']['H'] for m in range(1, 11)]), 100 * min(inpj[str(m)]['composition']['H'] for m in range(1, 11)), 100 * max(inpj[str(m)]['composition']['H'] for m in range(1, 11))))
