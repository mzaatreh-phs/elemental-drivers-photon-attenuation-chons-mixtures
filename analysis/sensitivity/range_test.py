#!/usr/bin/env python3
"""Range / distribution sensitivity of the elemental importance map.
Importance = Pratt share  s_i(E) = Delta_i * Cov(w_i, MAC) / Var(MAC)  with Delta_i = mu_i - mu_C  (sums to 1 for the exact linear mixture rule; handles correlated inputs).
Scenarios change the composition ensemble (ranges, distribution, N-S coupling); an analytic equal-variability case removes the range assumption altogether."""
import os, csv, json
import numpy as np
R = os.path.expanduser('~/CHON/recalc_2026-09-19'); OUT = R + '/analysis/sensitivity/'
rows = list(csv.DictReader(open(OUT + 'element_mac.csv'))); E = np.array([float(r['Energy_MeV']) for r in rows])
EL = ['C', 'H', 'N', 'O', 'S']; MU = np.array([[float(r[f'mu_{e}_cm2_per_g']) for r in rows] for e in EL])      # [5, nE]
inp = json.load(open(R + '/inputs/inputs_exact.json'))['mixtures']
W = np.array([[inp[str(m)]['composition'].get(e, 0.0) for e in EL] for m in range(1, 11)]); lo, hi = W.min(0), W.max(0); mid = (lo + hi) / 2; half = (hi - lo) / 2
rng = np.random.default_rng(20260919); NS = 20000
def box(scale, clip0=True, c_range=None, n=NS, tie_S=None):
    out = []
    while len(out) < n:
        x = rng.uniform(np.maximum(mid - scale * half, 0) if clip0 else mid - scale * half, mid + scale * half)          # x over [C,H,N,O,S]; C overwritten by closure
        if tie_S: x[4] = tie_S * x[2]
        c = 1 - x[1:].sum(); x[0] = c
        if (c_range is None and c > 0.2) or (c_range is not None and c_range[0] <= c <= c_range[1]): out.append(x)
    return np.array(out)
def dirich(kappa, n=NS):
    m = W.mean(0); return rng.dirichlet(m * kappa, size=n)
SC = {
 'A baseline (library box, x1)':        box(1.0, c_range=(lo[0], hi[0])),
 'B narrow (x0.5)':                     box(0.5),
 'C wide (x2)':                         box(2.0),
 'D N-S tied (BSA ratio S/N=0.114)':    box(1.0, tie_S=0.01883 / 0.16475),
 'E Dirichlet (kappa=300)':             dirich(300),
 'F broad organic space':               None,
}
# F: broad organic space (hydrocarbon-to-protein-to-sulfur-rich)
BL = np.array([0, 0.04, 0, 0, 0]); BH = np.array([1, 0.14, 0.25, 0.60, 0.05]); out = []
while len(out) < NS:
    x = rng.uniform(BL, BH); x[0] = 1 - x[1:].sum()
    if 0.2 <= x[0] <= 0.85: out.append(x)
SC['F broad organic space'] = np.array(out)
def pratt(X):
    mac = X @ MU; D = MU[1:] - MU[0]                                  # [4, nE]
    v = mac.var(0); cov = np.array([[np.cov(X[:, k], mac[:, j])[0, 1] for j in range(len(E))] for k in range(1, 5)])
    return D * cov / v                                                # [4, nE]  (H,N,O,S)
res = {k: pratt(v) for k, v in SC.items()}
# G: analytic equal-variability (independent inputs, equal SD): share ~ Delta^2
D = MU[1:] - MU[0]; res['G equal variability (analytic)'] = D ** 2 / (D ** 2).sum(0)
names = ['H', 'N', 'O', 'S']
with open(OUT + 'range_test_importance.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['scenario', 'Energy_MeV'] + [f'share_{n}' for n in names])
    for k, s in res.items():
        for j, e in enumerate(E): w.writerow([k, f'{e:g}'] + [f'{s[i, j]:.5f}' for i in range(4)])
def windows(s):
    dom = [names[int(np.argmax(s[:, j]))] for j in range(len(E))]; out = []; st = 0
    for j in range(1, len(E) + 1):
        if j == len(E) or dom[j] != dom[st]: out.append((dom[st], E[st], E[j - 1])); st = j
    return out
print('scenario ensembles: sd of (H,N,O,S) mass fractions, and corr(N,S):')
for k, X in SC.items(): print(f'  {k:38s} sd = {np.round(X[:, 1:].std(0), 4)}   corr(N,S)={np.corrcoef(X[:, 2], X[:, 4])[0, 1]:+.2f}')
print('\ndominant element by energy window (largest Pratt share):')
for k, s in res.items(): print(f'  {k:38s}', '; '.join(f'{d}: {a:g}-{b:g}' for d, a, b in windows(s)))
print('\nPratt share (%) of MAC variability at key energies  [H / N / O / S]')
KE = [0.01, 0.03, 0.05, 0.1, 1, 10, 15]
print('  ' + ' ' * 38 + ''.join(f'{e:>19g}' for e in KE))
for k, s in res.items():
    print(f'  {k:38s}' + ''.join('   ' + '/'.join(f'{100 * s[i, int(np.argmin(abs(E - e)))]:3.0f}' for i in range(4)) + '  ' for e in KE))
# robustness: fraction of the independent scenarios (A,B,C,E,F,G) agreeing on the dominant element at each energy
ind = [k for k in res if not k.startswith('D')]
agree = []
for j, e in enumerate(E):
    d = [names[int(np.argmax(res[k][:, j]))] for k in ind]; top = max(set(d), key=d.count); agree.append((e, top, d.count(top) / len(d)))
print('\nrobust dominant element (independent-input scenarios A,B,C,E,F,G): energy -> element (agreement)')
print('  ' + '; '.join(f'{e:g}:{t}({a:.0%})' for e, t, a in agree[::3]))
low = [a for e, t, a in agree if e <= 0.05]; mid_ = [a for e, t, a in agree if 0.1 <= e <= 3]; hi_ = [a for e, t, a in agree if e >= 10]
print(f'  agreement <=0.05 MeV: min {min(low):.0%}; 0.1-3 MeV: min {min(mid_):.0%}; >=10 MeV: min {min(hi_):.0%}')
