#!/usr/bin/env python3
"""Main analysis: pool seeds, check uncertainty, compare with XCOM, derive shielding parameters, spreads and correlations
with propagated noise, and write CSV / LaTeX / JSON outputs into analysis/g4/results/.
Works on partial data (only energies with all 10 mixtures x 2 seeds finished are used)."""
import os, sys, json, csv
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from g4analysis import *

os.makedirs(f'{OUT}/results', exist_ok=True)
RES = f'{OUT}/results/'
runs, problems = load_runs()
energies = complete_energies(runs)
if not energies:
    print('no complete energies yet'); sys.exit(0)
D = build_tables(runs, energies)
E = D['E']; nE = len(E)
xc = load_xcom_mac(); zeffd, zref = load_zeff()
NUM = {'n_jobs_loaded': len(runs), 'n_energies_complete': nE, 'problems': problems[:10]}

# ---- 0. raw-consistency: recomputed MAC must equal the code's own CSV value
dev = [abs(mac_from_counts(r['n_inc'], r['n_tr'], r['x'], r['rho'])[0] / r['csv_mac'] - 1) for r in runs.values()]
NUM['max_rel_diff_recomputed_vs_csv_MAC'] = float(max(dev))

# ---- 1. seed-to-seed check of the counting uncertainty
z = (D['mac1'] - D['mac2']) / np.sqrt(D['u1'] ** 2 + D['u2'] ** 2)
NUM['seed_check'] = {
    'n': int(z.size), 'mean_z': float(z.mean()), 'std_z': float(z.std(ddof=1)),
    'chi2': float((z ** 2).sum()), 'dof': int(z.size), 'p_value': float(1 - stats.chi2.cdf((z ** 2).sum(), z.size)),
    'frac_abs_z_gt2': float((np.abs(z) > 2).mean()),
    'by_band': {name: {'std_z': float(z[:, (E >= lo) & (E < hi)].std(ddof=1))}
                for name, lo, hi in (('<0.1 MeV', 0, 0.1), ('0.1-1 MeV', 0.1, 1), ('>=1 MeV', 1, 99)) if ((E >= lo) & (E < hi)).sum()}}
np.savetxt(RES + 'seed_check_z.csv', z, delimiter=',', header='rows=Mix1..10, cols=energies: ' + ','.join(f'{e:g}' for e in E))

# ---- 2. XCOM comparison
X = np.array([[xc[(m, round(e, 6))] for e in E] for m in MIXES])
dpct = 100 * (D['mac'] / X - 1)               # signed deviation, %
zx = (D['mac'] - X) / D['u']                  # significance in sigma
isval = np.array([round(e, 6) in [round(v, 6) for v in VALID] for e in E])
NUM['xcom'] = {
    'validation_max_abs_pct': float(np.abs(dpct[:, isval]).max()) if isval.any() else None,
    'validation_mean_abs_pct': float(np.abs(dpct[:, isval]).mean()) if isval.any() else None,
    'production_mean_abs_pct': float(np.abs(dpct[:, ~isval]).mean()) if (~isval).any() else None,
    'production_max_abs_pct': float(np.abs(dpct[:, ~isval]).max()) if (~isval).any() else None,
    'all_max_abs_pct': float(np.abs(dpct).max()),
    'frac_within_2sigma': float((np.abs(zx) < 2).mean()),
    'mean_signed_pct_by_energy': {f'{e:g}': float(dpct[:, j].mean()) for j, e in enumerate(E)},
    'max_pooled_rel_uncertainty_pct': float(100 * (D['u'] / D['mac']).max()),
    'min_pooled_rel_uncertainty_pct': float(100 * (D['u'] / D['mac']).min())}

with open(RES + 'mac_pooled.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['mixture', 'E_MeV', 'n_incident_pooled', 'n_transmitted_pooled', 'T', 'MAC_G4', 'u_MAC', 'MAC_XCOM',
                                   'dev_pct', 'dev_sigma', 'MAC_s1', 'MAC_s2', 'LAC', 'u_LAC', 'MFP_cm', 'HVL_cm', 'TVL_cm', 'RPE_1cm_pct', 'u_RPE_pct'])
    for i, m in enumerate(MIXES):
        for j, e in enumerate(E):
            lac = D['lac'][i, j]
            w.writerow([m, e, int(D['nP'][i, j]), int(D['kP'][i, j]), f"{D['T'][i, j]:.6f}", f"{D['mac'][i, j]:.6f}", f"{D['u'][i, j]:.6f}", f"{X[i, j]:.6f}",
                        f"{dpct[i, j]:+.3f}", f"{zx[i, j]:+.2f}", f"{D['mac1'][i, j]:.6f}", f"{D['mac2'][i, j]:.6f}", f"{lac:.6f}", f"{D['u_lac'][i, j]:.6f}",
                        f"{1 / lac:.5f}", f"{np.log(2) / lac:.5f}", f"{np.log(10) / lac:.5f}", f"{100 * (1 - D['T'][i, j]):.4f}", f"{100 * D['sigT'][i, j]:.4f}"])

# ---- 3. Table 3 (validation energies) as LaTeX in the manuscript layout
if isval.any():
    cols = [j for j, e in enumerate(E) if isval[j]]
    head = '& ' + ' & '.join(r'\multicolumn{3}{c}{%d keV}' % round(E[j] * 1000) for j in cols) + r' \\'
    cm = ''.join(r'\cmidrule(lr){%d-%d}' % (2 + 3 * k, 4 + 3 * k) for k in range(len(cols)))
    lines = [r'\begin{tabular}{c ' + ' '.join(['ccc'] * len(cols)) + '}', r'\toprule', head, cm,
             'Mix & ' + ' & '.join(['G4 & XCOM & $\\Delta\\%$'] * len(cols)) + r' \\', r'\midrule']
    for i, m in enumerate(MIXES):
        lines.append(f'{m} & ' + ' & '.join(f"{D['mac'][i, j]:.4f} & {X[i, j]:.4f} & {abs(dpct[i, j]):.3f}" for j in cols) + r' \\')
    lines += [r'\bottomrule', r'\end{tabular}']
    open(RES + 'table3_validation.tex', 'w').write('\n'.join(lines) + '\n')

# ---- 4. spreads with propagated noise (Table 4)
SEL = [0.01, 0.03, 0.1, 1.0, 10.0, 15.0]
selj = [int(np.argmin(np.abs(E - s))) for s in SEL if np.abs(E - s).min() < 1e-9]
t4 = []
for j in selj:
    row = {'E': float(E[j])}
    for q, mean, sd in (('MAC', D['mac'][:, j], D['u'][:, j]), ('LAC', D['lac'][:, j], D['u_lac'][:, j])):
        draws = mc_draw(mean, sd)
        sp = np.array([spread_pct(d) for d in draws]); lo, med, hi = interval(sp)
        p_max = float((draws.argmax(1) == mean.argmax()).mean()); p_min = float((draws.argmin(1) == mean.argmin()).mean())
        row[q] = dict(max_mix=int(np.argmax(mean)) + 1, min_mix=int(np.argmin(mean)) + 1, spread_pct=float(spread_pct(mean)),
                      mc_16=float(lo), mc_50=float(med), mc_84=float(hi), p_max_identity=p_max, p_min_identity=p_min,
                      xcom_spread_pct=float(spread_pct(X[:, j])) if q == 'MAC' else None)
    t4.append(row)
NUM['spreads'] = t4
json.dump(t4, open(RES + 'table4_spreads.json', 'w'), indent=1)
L = [r'\begin{tabular}{cccccc}', r'\toprule', r'Energy (MeV) & Max-MAC mix & Min-MAC mix & MAC spread (\%) & 68\% interval (\%) & LAC spread (\%) \\', r'\midrule']
for r in t4:
    L.append(f"{r['E']:g} & {r['MAC']['max_mix']} & {r['MAC']['min_mix']} & {r['MAC']['spread_pct']:.2f} & {r['MAC']['mc_16']:.2f}--{r['MAC']['mc_84']:.2f} & {r['LAC']['spread_pct']:.1f} \\\\")
L += [r'\bottomrule', r'\end{tabular}']
open(RES + 'table4_spread.tex', 'w').write('\n'.join(L) + '\n')

# ---- 5. correlations (Table 5): observed, noise-propagated interval, and noise-free XCOM theory
PRED = {'Oxygen': np.array([COMP[m]['O'] for m in MIXES]), 'Hydrogen': np.array([COMP[m]['H'] for m in MIXES]),
        'Carbon': np.array([COMP[m]['C'] for m in MIXES]), 'Nitrogen': np.array([COMP[m]['N'] for m in MIXES]),
        'Sulfur': np.array([COMP[m]['S'] for m in MIXES]), 'Density': np.array([RHO[m] for m in MIXES]),
        'BSA fraction': np.array([BSA[m] for m in MIXES])}
CORE = [0.01, 0.03, 0.05, 0.08, 0.1, 1.0, 10.0, 15.0]
corr = []
for pname, pv in PRED.items():
    for resp in ('MAC', 'LAC', 'Zeff'):
        for s in CORE:
            if np.abs(E - s).min() > 1e-9: continue
            j = int(np.argmin(np.abs(E - s)))
            if resp == 'MAC': obs = D['mac'][:, j]; sd = D['u'][:, j]; th = X[:, j]
            elif resp == 'LAC': obs = D['lac'][:, j]; sd = D['u_lac'][:, j]; th = X[:, j] * np.array([RHO[m] for m in MIXES])
            else: obs = np.array([zeffd[(m, round(E[j], 6))] for m in MIXES]); sd = None; th = obs
            r_obs = pearson(pv, obs)
            if sd is not None:
                rr = np.array([pearson(pv, d) for d in mc_draw(obs, sd, n=2000)]); lo, med, hi = interval(rr)
            else: lo = med = hi = r_obs
            corr.append(dict(pred=pname, resp=resp, E=float(E[j]), r_G4=r_obs, r16=float(lo), r50=float(med), r84=float(hi), r_xcom=pearson(pv, th)))
NUM['n_correlations'] = len(corr)
with open(RES + 'table5_correlations.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(corr[0].keys())); w.writeheader(); w.writerows(corr)
Es = [s for s in CORE if np.abs(E - s).min() < 1e-9]
L = [r'\begin{tabular}{l' + 'c' * len(Es) + '}', r'\toprule', 'Variable pair & ' + ' & '.join(f'{s:g}' for s in Es) + r' \\', r'\midrule']
for pname in PRED:
    for resp in ('MAC', 'LAC'):
        rr = {c['E']: c for c in corr if c['pred'] == pname and c['resp'] == resp}
        L.append(f'{pname} vs.\\ {resp} & ' + ' & '.join(f"{rr[s]['r_G4']:+.3f}" for s in Es if s in rr) + r' \\')
for pname in ('Oxygen', 'Sulfur'):
    rr = {c['E']: c for c in corr if c['pred'] == pname and c['resp'] == 'Zeff'}
    L.append(f'{pname} vs.\\ $Z_\\mathrm{{eff}}$ & ' + ' & '.join(f"{rr[s]['r_G4']:+.3f}" for s in Es if s in rr) + r' \\')
L += [r'\bottomrule', r'\end{tabular}']
open(RES + 'table5_correlation.tex', 'w').write('\n'.join(L) + '\n')

# ---- 5b. sign-change energies of O-MAC and S-MAC correlations over the whole complete grid (G4 and noise-free XCOM)
def sign_changes(resp_fn):
    rs = [pearson(PRED_, resp_fn(j)) for j in range(nE)]
    return [float(E[j]) for j in range(1, nE) if rs[j - 1] * rs[j] < 0], rs
sc = {}
for pname in ('Oxygen', 'Sulfur', 'Hydrogen'):
    PRED_ = PRED[pname]
    g4, r_g4 = sign_changes(lambda j: D['mac'][:, j]); th, r_th = sign_changes(lambda j: X[:, j])
    sc[pname] = dict(G4_sign_change_between=g4, XCOM_sign_change_between=th)
NUM['sign_changes_MAC_correlation'] = sc

# ---- 6. derived parameters (ranges for text) at selected energies
rng = {}
for s in SEL:
    if np.abs(E - s).min() > 1e-9: continue
    j = int(np.argmin(np.abs(E - s)))
    rng[f'{s:g}'] = dict(MAC=[float(D['mac'][:, j].min()), float(D['mac'][:, j].max())], LAC=[float(D['lac'][:, j].min()), float(D['lac'][:, j].max())],
                         HVL_cm=[float((np.log(2) / D['lac'][:, j]).min()), float((np.log(2) / D['lac'][:, j]).max())],
                         RPE_1cm_pct=[float(100 * (1 - D['T'][:, j]).min()), float(100 * (1 - D['T'][:, j]).max())],
                         thickness_saving_pct_maxLAC_vs_minLAC=float(100 * (1 - D['lac'][:, j].min() / D['lac'][:, j].max())),
                         areal_mass_ratio_maxLACmix_vs_minLACmix=float((1 / D['mac'][int(np.argmax(D['lac'][:, j])), j]) / (1 / D['mac'][int(np.argmin(D['lac'][:, j])), j])))
NUM['ranges_at_selected_energies'] = rng

json.dump(NUM, open(RES + 'numbers_for_text.json', 'w'), indent=1, default=float)
print(f"complete energies: {nE} of {len(GRID)}   jobs loaded: {len(runs)}   problems: {len(problems)}")
print('recomputed-vs-csv MAC max rel diff: %.2e' % NUM['max_rel_diff_recomputed_vs_csv_MAC'])
sc_ = NUM['seed_check']; print('seed check: std(z)=%.3f (expect ~1), mean z=%+.3f, chi2/dof=%.3f, p=%.3f, |z|>2: %.1f%%' % (sc_['std_z'], sc_['mean_z'], sc_['chi2'] / sc_['dof'], sc_['p_value'], 100 * sc_['frac_abs_z_gt2']))
xr = NUM['xcom']; print('XCOM: mean |dev| prod %.3f%%, max |dev| %.3f%%, within 2 sigma %.0f%%, pooled MAC rel. unc. %.2f-%.2f%%' % (xr['production_mean_abs_pct'] or float('nan'), xr['all_max_abs_pct'], 100 * xr['frac_within_2sigma'], xr['min_pooled_rel_uncertainty_pct'], xr['max_pooled_rel_uncertainty_pct']))
