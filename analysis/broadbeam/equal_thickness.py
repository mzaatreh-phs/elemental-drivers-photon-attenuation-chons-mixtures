#!/usr/bin/env python3
"""Item 2: total (uncollided + scattered) exposure response  T_broad = B(E,d) exp(-d)  at EQUAL PHYSICAL THICKNESS and at EQUAL AREAL MASS.
G-P buildup (ANSI/ANS-6.4.3 via epixs_cli, mixture Zeq path) with XCOM total MAC (incl. coherent); d = mu t = MAC*rho*t (thickness) or MAC*m (areal mass).
G-P is defined for 0.5 <= d <= 40 mfp; outside that range the entry is left empty.  Point isotropic source, infinite homogeneous medium, exposure response."""
import sys, os, csv, json
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.expanduser('~/epixs_cli')); sys.path.insert(0, os.path.expanduser('~/epixs_cli/epixs'))
from zeff import zeq
from gp_calculator import material_gp_coeffs, gp_buildup_factor
R = os.path.expanduser('~/CHON/recalc_2026-09-19'); OUT = R + '/analysis/broadbeam/'
inp = json.load(open(R + '/inputs/inputs_exact.json')); mix = inp['mixtures']
E = [e for e in sorted(inp['production_energies'] + inp['validation_energies']) if e >= 0.015]
MIXES = list(range(1, 11)); rho = {m: mix[str(m)]['density'] for m in MIXES}
xc = {}
for r in csv.DictReader(open(R + '/analysis/xcom/process_fractions.csv')):
    if r['set'] == 'new': xc[(int(r['material'][3:]), round(float(r['E_MeV']), 6))] = float(r['mac_total_cm2_g'])
g4 = {}; g4u = {}
for r in csv.DictReader(open(R + '/analysis/g4/results/mac_pooled.csv')):
    g4[(int(r['mixture']), round(float(r['E_MeV']), 6))] = float(r['MAC_G4']); g4u[(int(r['mixture']), round(float(r['E_MeV']), 6))] = float(r['u_MAC']) / float(r['MAC_G4'])
COEF = {}
for m in MIXES:
    comp = {k: v for k, v in mix[str(m)]['composition'].items() if v > 0}
    for e in E:
        z, z1, z2 = zeq(comp, e); COEF[(m, e)] = material_gp_coeffs(e, z, z1, z2, 'exposure')
def B(m, e, d):
    b, c, a, Xk, dd = COEF[(m, e)]; return gp_buildup_factor(d, b, c, a, Xk, dd)
THK = [1, 2, 5, 10, 20, 50, 100]          # cm
AM = [1, 2, 5, 10, 20, 50, 100]           # g/cm2
def valid(d): return 0.5 <= d <= 40
rows = []
for kind, vals in (('thickness_cm', THK), ('areal_mass_g_cm2', AM)):
    for v in vals:
        for e in E:
            for m in MIXES:
                mac = xc[(m, round(e, 6))]; macg = g4[(m, round(e, 6))]
                d = mac * rho[m] * v if kind == 'thickness_cm' else mac * v
                dg = macg * rho[m] * v if kind == 'thickness_cm' else macg * v
                if valid(d):
                    b = B(m, e, d); rows.append([kind, v, e, m, d, b, np.exp(-d), b * np.exp(-d), (B(m, e, dg) * np.exp(-dg)) if valid(dg) else ''])
                else: rows.append([kind, v, e, m, d, '', '', '', ''])
with open(OUT + 'broadbeam_long.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['condition', 'value', 'E_MeV', 'mixture', 'd_mfp', 'B_exposure', 'T_primary', 'T_broad', 'T_broad_using_Geant4_MAC']); w.writerows(rows)
# ---- summary per (condition, value, energy): spreads, best/worst, rank agreement
summ = []
for kind, vals in (('thickness_cm', THK), ('areal_mass_g_cm2', AM)):
    for v in vals:
        for e in E:
            sel = [r for r in rows if r[0] == kind and r[1] == v and r[2] == e and r[5] != '']
            if len(sel) < 10: continue
            Tp = np.array([r[6] for r in sel]); Tb = np.array([r[7] for r in sel]); mx = [r[3] for r in sel]; d = np.array([r[4] for r in sel])
            sp_p = 100 * (Tp.max() / Tp.min() - 1); sp_b = 100 * (Tb.max() / Tb.min() - 1)
            rho_s = stats.spearmanr(Tp, Tb).statistic
            summ.append([kind, v, e, d.min(), d.max(), sp_p, sp_b, mx[int(np.argmin(Tp))], mx[int(np.argmin(Tb))], mx[int(np.argmax(Tb))], rho_s])
with open(OUT + 'broadbeam_summary.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['condition', 'value', 'E_MeV', 'd_min_mfp', 'd_max_mfp', 'spread_primary_pct', 'spread_broad_pct', 'best_mix_primary', 'best_mix_broad', 'worst_mix_broad', 'spearman_primary_vs_broad'])
    w.writerows([[s[0], s[1], f'{s[2]:g}'] + [f'{x:.5g}' if isinstance(x, float) else x for x in s[3:]] for s in summ])
# ---- wide tables for plotting: key conditions
def wide(name, kind, v, col):
    hdr = ['Energy_MeV'] + [f'Mix{m}' for m in MIXES]; out = [hdr]
    for e in E:
        r = {rr[3]: rr[col] for rr in rows if rr[0] == kind and rr[1] == v and rr[2] == e}
        out.append([f'{e:g}'] + [f'{r[m]:.6g}' if r[m] != '' else '' for m in MIXES])
    csv.writer(open(OUT + name + '.csv', 'w', newline='')).writerows(out)
for kind, v, tag in (('thickness_cm', 10, 't10cm'), ('thickness_cm', 20, 't20cm'), ('areal_mass_g_cm2', 20, 'm20gcm2'), ('areal_mass_g_cm2', 50, 'm50gcm2')):
    wide(f'wide_T_broad_{tag}', kind, v, 7); wide(f'wide_T_primary_{tag}', kind, v, 6); wide(f'wide_B_{tag}', kind, v, 5)
# ---- print key findings
def show(kind, v):
    S = [s for s in summ if s[0] == kind and s[1] == v]
    print(f'\n{kind} = {v}: energies with valid G-P range: {len(S)}')
    print('  E(MeV)  d range(mfp)   spread prim%  spread broad%  best(prim) best(broad) worst(broad)  Spearman')
    for s in S[::max(1, len(S) // 12)]:
        print(f'  {s[2]:6g}  {s[3]:5.1f}-{s[4]:5.1f}     {s[5]:9.2f}   {s[6]:9.2f}       {s[7]:>3}       {s[8]:>3}        {s[9]:>3}       {s[10]:+.2f}')
    flips = [s for s in S if s[7] != s[8]]
    print(f'  best mixture differs between primary-only and with buildup at {len(flips)} of {len(S)} energies', ('(' + ', '.join(f'{s[2]:g}' for s in flips[:10]) + ')') if flips else '')
for kind, v in (('thickness_cm', 10), ('areal_mass_g_cm2', 20)): show(kind, v)
