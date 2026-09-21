#!/usr/bin/env python3
"""Item 1: elemental sensitivity map + decorrelated test ensemble (XCOM mixture rule, exact CHONS compositions).
MAC(E) = sum_i w_i mu_i(E) (mixture rule). With closure sum w = 1 and carbon as reference:
    MAC = mu_C + sum_{i in H,N,O,S} w_i (mu_i - mu_C)   -> contrast  Delta_i(E) = mu_i(E) - mu_C(E)."""
import sys, os, csv, json
import numpy as np
sys.path.insert(0, os.path.expanduser('~/epixs_cli'))
from epixs.xcom_engine import mass_attenuation
R = os.path.expanduser('~/CHON/recalc_2026-09-19'); OUT = R + '/analysis/sensitivity/'
inp = json.load(open(R + '/inputs/inputs_exact.json')); mix = inp['mixtures']
E = np.array(sorted(inp['production_energies'] + inp['validation_energies']))
EL = ['C', 'H', 'N', 'O', 'S']
fine = np.geomspace(0.005, 15.0, 700)
def mu(el, en): return np.array(mass_attenuation(el, list(en))['total_with_coherent'])
MU = {e: mu(e, E) for e in EL}; MUf = {e: mu(e, fine) for e in EL}
# ---- library compositions (ten mixtures)
W = np.array([[mix[str(m)]['composition'].get(e, 0.0) for e in EL] for m in range(1, 11)])     # columns C,H,N,O,S
mac_lib = W @ np.array([MU[e] for e in EL])                                                        # [10, nE]
# ---- element contrasts vs carbon and zero crossings (exchange crossovers)
def crossings(x, y):
    out = []
    for i in range(len(y) - 1):
        if y[i] == 0: out.append(x[i])
        elif y[i] * y[i + 1] < 0: out.append(float(np.exp(np.log(x[i]) + (0 - y[i]) / (y[i + 1] - y[i]) * (np.log(x[i + 1]) - np.log(x[i])))))
    return out
cross = {i: crossings(fine, MUf[i] - MUf['C']) for i in ('H', 'N', 'O', 'S')}
# ---- mean-library contributions and per-1wt% exchange sensitivity
wbar = W.mean(0); macbar = wbar @ np.array([MU[e] for e in EL])
contrib = {e: wbar[k] * MU[e] / macbar * 100 for k, e in enumerate(EL)}          # % of MAC carried by each element
sens = {i: (MU[i] - MU['C']) / macbar * 1.0 for i in ('H', 'N', 'O', 'S')}      # % MAC change per +1 wt% i (replacing C)
with open(OUT + 'element_mac.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['Energy_MeV'] + [f'mu_{e}_cm2_per_g' for e in EL])
    for j, en in enumerate(E): w.writerow([f'{en:g}'] + [f'{MU[e][j]:.8g}' for e in EL])
with open(OUT + 'exchange_sensitivity_and_contribution.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['Energy_MeV'] + [f'pct_MAC_change_per_1wt%_{i}_replacing_C' for i in ('H', 'N', 'O', 'S')] + [f'contribution_{e}_pct_of_MAC_library_mean' for e in EL])
    for j, en in enumerate(E): w.writerow([f'{en:g}'] + [f'{sens[i][j]:.6g}' for i in ('H', 'N', 'O', 'S')] + [f'{contrib[e][j]:.6g}' for e in EL])
# ---- decorrelated ensemble inside the library's composition ranges
rng = np.random.default_rng(20260919); lo = W.min(0); hi = W.max(0)
ens = []
while len(ens) < 4000:
    h = rng.uniform(lo[1], hi[1]); n = rng.uniform(lo[2], hi[2]); o = rng.uniform(lo[3], hi[3]); s = rng.uniform(lo[4], hi[4]); c = 1 - h - n - o - s
    if lo[0] <= c <= hi[0]: ens.append([c, h, n, o, s])
ens = np.array(ens)
def corrmat(X): return np.corrcoef(X[:, 1:].T)          # H,N,O,S
def vif(X):
    Z = X[:, 1:]; out = []
    for k in range(Z.shape[1]):
        y = Z[:, k]; A = np.delete(Z, k, 1); A = np.c_[np.ones(len(A)), A]
        beta, *_ = np.linalg.lstsq(A, y, rcond=None); r2 = 1 - ((y - A @ beta) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        out.append(1 / (1 - r2) if r2 < 1 - 1e-12 else float('inf'))
    return out
stats = {'library_corr_HNOS': corrmat(W).round(3).tolist(), 'ensemble_corr_HNOS': corrmat(ens).round(3).tolist(),
         'library_VIF_HNOS': [float(x) if np.isfinite(x) else 'inf' for x in vif(W)], 'ensemble_VIF_HNOS': [round(x, 2) for x in vif(ens)],
         'library_corr_N_S': float(np.corrcoef(W[:, 2], W[:, 4])[0, 1]), 'ensemble_corr_N_S': float(np.corrcoef(ens[:, 2], ens[:, 4])[0, 1]),
         'library_design_rank_of_[H,N,O,S]': int(np.linalg.matrix_rank(W[:, 1:] - W[:, 1:].mean(0))), 'ensemble_design_rank': int(np.linalg.matrix_rank(ens[:, 1:] - ens[:, 1:].mean(0)))}
MUmat = np.array([MU[e] for e in EL])                       # [5, nE]
mac_ens = ens @ MUmat                                        # [n, nE]
def rcol(X, M):  # correlation of each predictor (H,N,O,S) with MAC at each energy
    return np.array([[np.corrcoef(X[:, k], M[:, j])[0, 1] if X[:, k].std() > 0 else np.nan for j in range(M.shape[1])] for k in range(1, 5)])
r_lib = rcol(W, mac_lib); r_ens = rcol(ens, mac_ens)
# standardized (importance) coefficients on the ensemble: beta_i = Delta_i * sd(w_i) / sd(MAC)
beta = np.array([(MU[i] - MU['C']) * ens[:, k].std() / mac_ens.std(0) for k, i in zip(range(1, 5), ('H', 'N', 'O', 'S'))])
imp = np.abs(beta) / np.abs(beta).sum(0) * 100
with open(OUT + 'ensemble_vs_library_correlation_and_importance.csv', 'w', newline='') as f:
    w = csv.writer(f); hdr = ['Energy_MeV']
    for i in ('H', 'N', 'O', 'S'): hdr += [f'r_library_{i}', f'r_ensemble_{i}', f'beta_std_{i}', f'importance_pct_{i}']
    w.writerow(hdr)
    for j, en in enumerate(E):
        row = [f'{en:g}']
        for k in range(4): row += [f'{r_lib[k, j]:.4f}', f'{r_ens[k, j]:.4f}', f'{beta[k, j]:.4f}', f'{imp[k, j]:.2f}']
        w.writerow(row)
json.dump({'exchange_zero_crossings_MeV': cross, 'predictor_statistics': stats}, open(OUT + 'summary.json', 'w'), indent=1)
# dominant element by energy (importance) -> windows
names = ['H', 'N', 'O', 'S']; dom = [names[int(np.argmax(imp[:, j]))] for j in range(len(E))]
win = []; start = 0
for j in range(1, len(E) + 1):
    if j == len(E) or dom[j] != dom[start]: win.append((dom[start], float(E[start]), float(E[j - 1]))); start = j
json.dump(win, open(OUT + 'dominant_element_windows.json', 'w'))
print('exchange zero-crossings (MeV) where replacing C by X leaves MAC unchanged:'); [print(f'  {i}: ', ', '.join(f'{x*1000:.1f} keV' if x < 1 else f'{x:.2f} MeV' for x in c) or 'none') for i, c in cross.items()]
print('\nlibrary vs decorrelated ensemble: corr(N,S) library %.3f, ensemble %.3f ; design rank library %d, ensemble %d ; VIF library %s ; ensemble %s' % (stats['library_corr_N_S'], stats['ensemble_corr_N_S'], stats['library_design_rank_of_[H,N,O,S]'], stats['ensemble_design_rank'], stats['library_VIF_HNOS'], stats['ensemble_VIF_HNOS']))
print('\ndominant element (importance) by energy window:'); [print(f'  {d}: {a:g}-{b:g} MeV') for d, a, b in win]
print('\n E(MeV)  sensitivity % MAC per +1 wt% (replacing C):  H     N     O     S   | importance % H/N/O/S | r_lib O  r_ens O | r_lib H  r_ens H')
for en in (0.01, 0.02, 0.03, 0.05, 0.1, 0.3, 1, 3, 10, 15):
    j = int(np.argmin(abs(E - en)))
    print(f' {en:6g}                                          {sens["H"][j]:+6.2f} {sens["N"][j]:+6.2f} {sens["O"][j]:+6.2f} {sens["S"][j]:+6.2f} | {imp[0,j]:4.0f} {imp[1,j]:4.0f} {imp[2,j]:4.0f} {imp[3,j]:4.0f}   | {r_lib[2,j]:+.2f} {r_ens[2,j]:+.2f}  | {r_lib[0,j]:+.2f} {r_ens[0,j]:+.2f}')
