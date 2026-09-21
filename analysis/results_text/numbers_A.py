#!/usr/bin/env python3
"""Numbers for Results 3.2-3.5 (composition sensitivity, MAC/LAC, MFP/HVL/TVL, RPE)."""
import os, sys, csv, json
import numpy as np
sys.path.insert(0, os.path.expanduser('~/CHON/recalc_2026-09-19/analysis/g4'))
from g4analysis import *
runs, _ = load_runs(); energies = complete_energies(runs); D = build_tables(runs, energies); E = np.array(D['E'])
xc = load_xcom_mac(); X = np.array([[xc[(m, round(e, 6))] for e in E] for m in MIXES])
mac, lac, u, uT = D['mac'], D['lac'], D['u'], D['T']
rho = np.array([RHO[m] for m in MIXES]); W = {el: np.array([COMP[m].get(el, 0.0) for m in MIXES]) for el in 'CHONS'}
ix = lambda e: int(np.argmin(abs(E - e)))
def rng(a, fmt='%.4g'): return fmt % a.min() + '--' + fmt % a.max()
print('== library correlations among elements & density')
for a, b in (('O', 'H'), ('O', 'C'), ('O', 'N'), ('O', 'S'), ('H', 'C'), ('C', 'N')): print(' r(%s,%s)=%+.3f' % (a, b, np.corrcoef(W[a], W[b])[0, 1]), end=';')
print(' r(O,rho)=%+.3f r(H,rho)=%+.3f r(C,rho)=%+.3f r(N,rho)=%+.3f' % tuple(np.corrcoef(W[k], rho)[0, 1] for k in 'OHCN'))
print('== MAC (Geant4) min/max per selected energy: value range, mixtures, spread%; LAC; HVL; RPE')
for e in (0.01, 0.03, 0.1, 1, 10, 15):
    j = ix(e); m = mac[:, j]; l = lac[:, j]; r = 100 * (1 - uT[:, j])
    print(' E=%g MAC %s (max Mix%d min Mix%d) spread %.2f%% | LAC %s (max Mix%d min Mix%d) spread %.1f%% | HVL %s | MFP %s | TVL %s | RPE %s (max Mix%d min Mix%d) | u_MAC %.2f-%.2f%%' % (
        e, rng(m, '%.5g'), MIXES[m.argmax()], MIXES[m.argmin()], spread_pct(m), rng(l, '%.5g'), MIXES[l.argmax()], MIXES[l.argmin()], spread_pct(l),
        rng(np.log(2) / l), rng(1 / l), rng(np.log(10) / l), rng(r, '%.2f'), MIXES[r.argmax()], MIXES[r.argmin()], 100 * (u[:, j] / mac[:, j]).min(), 100 * (u[:, j] / mac[:, j]).max()))
print('== LAC ordering over the whole grid: densest/least dense mixture = largest/smallest LAC at how many of', len(E), 'energies')
print(' largest LAC is Mix10 at', int(sum(lac[:, j].argmax() == 9 for j in range(len(E)))), ' smallest LAC is Mix2 at', int(sum(lac[:, j].argmin() == 1 for j in range(len(E)))))
sp = np.array([spread_pct(lac[:, j]) for j in range(len(E))]); print(' LAC spread over grid: %.1f - %.1f%% (min at E=%g, max at E=%g)' % (sp.min(), sp.max(), E[sp.argmin()], E[sp.argmax()]))
sm = np.array([spread_pct(mac[:, j]) for j in range(len(E))]); print(' MAC spread over grid: %.2f - %.2f%%; XCOM MAC spread %.2f - %.2f%%' % (sm.min(), sm.max(), min(spread_pct(X[:, j]) for j in range(len(E))), max(spread_pct(X[:, j]) for j in range(len(E)))))
for e in (0.06, 0.08, 0.1, 0.2, 0.5, 1, 3):
    j = ix(e); print('  spread MAC G4 %.2f XCOM %.2f at %g MeV' % (spread_pct(mac[:, j]), spread_pct(X[:, j]), e))
rl = np.array([pearson(list(rho), list(lac[:, j])) for j in range(len(E))]); print(' r(rho,LAC) over grid %.3f - %.3f ; at 6 MeV %.3f ; at 15 MeV %.3f' % (rl.min(), rl.max(), rl[ix(6)], rl[ix(15)]))
for nm, f in (('MFP', lambda j: 1 / lac[:, j]), ('RPE', lambda j: 100 * (1 - uT[:, j]))):
    r = np.array([pearson(list(rho), list(f(j))) for j in range(len(E))]); print(' r(rho,%s) over grid %.3f to %.3f; at 0.03,0.1,1,10,15: %s' % (nm, r.min(), r.max(), [round(r[ix(e)], 3) for e in (0.03, 0.1, 1, 10, 15)]))
print('== equal-transmission comparison densest(Mix10) vs least dense (Mix2)')
for e in (0.03, 0.1, 1):
    j = ix(e); l10, l2 = lac[9, j], lac[1, j]; m10, m2 = mac[9, j], mac[1, j]
    ul10, ul2 = D['u_lac'][9, j], D['u_lac'][1, j]; um10, um2 = u[9, j], u[1, j]
    th = 100 * (1 - l2 / l10); am = 100 * (m2 / m10 - 1)   # thickness Mix10 needs vs Mix2 (ratio l2/l10); areal mass Mix10 needs vs Mix2 = m2/m10
    uam = 100 * (m2 / m10) * np.hypot(um2 / m2, um10 / m10)
    print(' E=%g: Mix10 needs %.1f%% less thickness than Mix2; areal mass of Mix10 is %+.2f%% (+-%.2f) relative to Mix2; MAC10=%.4f MAC2=%.4f' % (e, th, am, uam, m10, m2))
print(' H fraction Mix2 %.4f, Mix10 %.4f' % (W['H'][1], W['H'][9]))
print('== MAC ordering at 0.1 MeV (sorted mixtures by MAC):', [MIXES[i] for i in np.argsort(-mac[:, ix(0.1)])])
print('== XCOM MAC ordering at 0.1 MeV:', [MIXES[i] for i in np.argsort(-X[:, ix(0.1)])])
print('== HVL over the grid: range min-max at first/last:', rng(np.log(2) / lac[:, 0]), rng(np.log(2) / lac[:, -1]))
print('== RPE: highest RPE mixture at each of the grid energies == Mix10?', int(sum((uT[:, j].argmin() == 9) for j in range(len(E)))), 'of', len(E), '; lowest == Mix2', int(sum((uT[:, j].argmax() == 1) for j in range(len(E)))))
