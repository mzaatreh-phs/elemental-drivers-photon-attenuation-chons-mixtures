#!/usr/bin/env python3
"""Compare Geant4 controlled-exchange results with the XCOM mixture-rule prediction (relative MAC change vs the base composition)."""
import sys, os, csv, glob, json
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.expanduser('~/epixs_cli'))
from epixs.xcom_engine import mass_attenuation
R = os.path.expanduser('~/CHON/recalc_2026-09-19'); OUT = R + '/analysis/sensitivity/'
des = json.load(open(R + '/inputs/exchange_design.json')); rho = des['density_g_cm3']; EN = des['energies_MeV']; variants = list(des['design'])
raw = {}
for d in glob.glob(R + '/runs/exchange/s[12]_*_E*'):
    if not os.path.exists(d + '/done.flag'): continue
    name = os.path.basename(d); seed, var, ei = name.split('_')[0], name.split('_')[1], int(name.split('_E')[1])
    r = list(csv.DictReader(open(glob.glob(d + '/MAC_results_*.csv')[0])))[0]
    raw[(seed, var, ei)] = (int(r['n_incident']), int(r['n_transmitted_strict']))
def mac(nk):  # pooled over both seeds
    n = sum(v[0] for v in nk); k = sum(v[1] for v in nk); T = k / n
    return -np.log(T) / rho, np.sqrt((1 - T) / (n * T)) / rho
pred = {v: np.array(mass_attenuation({e: f for e, f in des['design'][v]['composition'].items() if f > 0}, EN)['total_with_coherent']) for v in variants}
rows = []; zs = []
for i, E in enumerate(EN):
    have = all((s, v, i) in raw for s in ('s1', 's2') for v in variants)
    if not have: continue
    G = {v: mac([raw[('s1', v, i)], raw[('s2', v, i)]]) for v in variants}
    for v in variants[1:]:
        d_g = G[v][0] / G['base'][0] - 1; u = np.hypot(G[v][1] / G[v][0], G['base'][1] / G['base'][0]) * (1 + d_g)
        d_x = pred[v][i] / pred['base'][i] - 1; z = (d_g - d_x) / u
        rows.append([E, v, 100 * d_g, 100 * u, 100 * d_x, z, G['base'][0], G[v][0], pred['base'][i], pred[v][i]]); zs.append(z)
    # additivity test: dALL vs sum of single exchanges (linear mixture rule => equal)
    s_g = sum(G[v][0] / G['base'][0] - 1 for v in ('dO', 'dN', 'dS', 'dH')); a_g = G['dALL'][0] / G['base'][0] - 1
    rows.append([E, 'additivity(sum singles vs ALL)', 100 * a_g, 0, 100 * s_g, np.nan, np.nan, np.nan, np.nan, np.nan])
with open(OUT + 'exchange_confirmation.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['E_MeV', 'variant', 'dMAC_Geant4_pct', 'u_pct', 'dMAC_XCOM_mixture_rule_pct', 'z', 'MAC_base_G4', 'MAC_variant_G4', 'MAC_base_XCOM', 'MAC_variant_XCOM'])
    for r in rows: w.writerow([f'{r[0]:g}', r[1]] + [f'{x:.5g}' for x in r[2:]])
zs = np.array(zs)
print(f'{len(raw)} of {len(EN)*2*len(variants)} jobs finished; {len(zs)} (variant, energy) comparisons')
if len(zs):
    print('agreement of Geant4 exchange effect with mixture rule: mean z %+.2f, std z %.2f, chi2/dof %.2f, |z|>2: %d of %d' % (zs.mean(), zs.std(ddof=1), (zs ** 2).sum() / len(zs), (abs(zs) > 2).sum(), len(zs)))
    print('\n E(MeV) variant   Geant4 dMAC%   +-     XCOM dMAC%    z')
    for r in rows:
        if r[1].startswith('additivity'): print(f' {r[0]:6g} ADDITIVITY: ALL {r[2]:+.2f}%  vs sum of singles {r[4]:+.2f}%')
        else: print(f' {r[0]:6g} {r[1]:6s}   {r[2]:+9.2f}  {r[3]:5.2f}   {r[4]:+9.2f}   {r[5]:+.1f}')
