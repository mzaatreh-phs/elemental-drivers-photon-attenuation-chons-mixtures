#!/usr/bin/env python3
"""Residual correlations, standard deviations and closure conditions of the seven range-test scenarios (A-G) of range_test.py, regenerated with the identical seed and
sampling code, plus a check that the Pratt variance shares reproduce range_test_importance.csv."""
import os, csv, json
import numpy as np
R = os.path.expanduser('~/CHON/recalc_2026-09-19'); OUT = R + '/analysis/sensitivity/'
rows = list(csv.DictReader(open(OUT + 'element_mac.csv'))); E = np.array([float(r['Energy_MeV']) for r in rows])
EL = ['C', 'H', 'N', 'O', 'S']; MU = np.array([[float(r[f'mu_{e}_cm2_per_g']) for r in rows] for e in EL])
inp = json.load(open(R + '/inputs/inputs_exact.json'))['mixtures']
W = np.array([[inp[str(m)]['composition'].get(e, 0.0) for e in EL] for m in range(1, 11)]); lo, hi = W.min(0), W.max(0); mid = (lo + hi) / 2; half = (hi - lo) / 2
rng = np.random.default_rng(20260919); NS = 20000
def box(scale, clip0=True, c_range=None, n=NS, tie_S=None):
    out = []
    while len(out) < n:
        x = rng.uniform(np.maximum(mid - scale * half, 0) if clip0 else mid - scale * half, mid + scale * half)
        if tie_S: x[4] = tie_S * x[2]
        c = 1 - x[1:].sum(); x[0] = c
        if (c_range is None and c > 0.2) or (c_range is not None and c_range[0] <= c <= c_range[1]): out.append(x)
    return np.array(out)
def dirich(kappa, n=NS): m = W.mean(0); return rng.dirichlet(m * kappa, size=n)
SC = {'A baseline (library box, x1)': box(1.0, c_range=(lo[0], hi[0])), 'B narrow (x0.5)': box(0.5), 'C wide (x2)': box(2.0),
      'D N-S tied (BSA ratio S/N=0.114)': box(1.0, tie_S=0.01883 / 0.16475), 'E Dirichlet (kappa=300)': dirich(300), 'F broad organic space': None}
BL = np.array([0, 0.04, 0, 0, 0]); BH = np.array([1, 0.14, 0.25, 0.60, 0.05]); out = []
while len(out) < NS:
    x = rng.uniform(BL, BH); x[0] = 1 - x[1:].sum()
    if 0.2 <= x[0] <= 0.85: out.append(x)
SC['F broad organic space'] = np.array(out)
def pratt(X):
    mac = X @ MU; D = MU[1:] - MU[0]; v = mac.var(0)
    cov = np.array([[np.cov(X[:, k], mac[:, j])[0, 1] for j in range(len(E))] for k in range(1, 5)]); return D * cov / v
# sanity: reproduce the stored shares
stored = list(csv.DictReader(open(OUT + 'range_test_importance.csv'))); mx = 0
for name, X in SC.items():
    sh = pratt(X)
    for j, e in enumerate(E):
        r = [x for x in stored if x['scenario'] == name and abs(float(x['Energy_MeV']) - e) < 1e-9][0]
        mx = max(mx, max(abs(sh[k, j] - float(r['share_' + el])) for k, el in enumerate('HNOS')))
print('max |regenerated Pratt share - stored share| over all scenarios/energies: %.2e' % mx)
names = ['H', 'N', 'O', 'S']; pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
with open(OUT + 'scenario_correlations.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['scenario', 'n', 'carbon_min', 'carbon_max'] + [f'sd_{n}' for n in names] + ['r_' + names[a] + names[b] for a, b in pairs] + ['max_abs_offdiag_r'])
    for name, X in SC.items():
        C = np.corrcoef(X[:, 1:].T); r = [C[a, b] for a, b in pairs]; sd = X[:, 1:].std(0)
        w.writerow([name, len(X), '%.4f' % X[:, 0].min(), '%.4f' % X[:, 0].max()] + ['%.4f' % v for v in sd] + ['%.3f' % v for v in r] + ['%.3f' % np.nanmax(np.abs(r))])
        print('%-36s n=%d C in [%.3f,%.3f]  sd(H,N,O,S)=%s  r: %s' % (name, len(X), X[:, 0].min(), X[:, 0].max(), np.round(sd, 4), {names[a] + names[b]: round(C[a, b], 2) for a, b in pairs}))
