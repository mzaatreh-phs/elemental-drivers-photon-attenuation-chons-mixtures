#!/usr/bin/env python3
"""Moisture sensitivity: each dry mixture is mixed with 10 wt% water (wet-mass basis: 0.9 dry + 0.1 H2O), as reported for the conditioned precursor powders
of Campbell et al. 2025.  XCOM (mixture rule, via EpiXS-CLI), Zeff (direct), crossover, and G-P exposure buildup are recomputed for the wet compositions.
Density: (a) additive volume with water at 1.00 g/cm3, (b) unchanged."""
import sys, os, csv, json, math
import numpy as np
sys.path.insert(0, os.path.expanduser('~/epixs_cli')); sys.path.insert(0, os.path.expanduser('~/epixs_cli/epixs'))
from zeff import zeff, zeq
from gp_calculator import material_gp_coeffs, gp_buildup_factor
from epixs.xcom_engine import mass_attenuation
R = os.path.expanduser('~/CHON/recalc_2026-09-19'); O = R + '/analysis/moisture/'
inp = json.load(open(R + '/inputs/inputs_exact.json')); mix = inp['mixtures']; E38 = sorted(inp['production_energies'] + inp['validation_energies'])
WATER = {'H': 0.111894, 'O': 0.888106}; FW = 0.10; RHOW = 1.0
M = list(range(1, 11)); dry = {m: mix[str(m)]['composition'] for m in M}; rho = {m: mix[str(m)]['density'] for m in M}
def wetc(c):
    out = {e: (1 - FW) * c.get(e, 0.0) + FW * WATER.get(e, 0.0) for e in set(c) | set(WATER)}; return {k: v for k, v in out.items() if v > 0}
wet = {m: wetc(dry[m]) for m in M}; rho_wet = {m: 1 / ((1 - FW) / rho[m] + FW / RHOW) for m in M}
def mac(c): return np.array(mass_attenuation(c, E38)['total_with_coherent'])
MD = np.array([mac(dry[m]) for m in M]); MW = np.array([mac(wet[m]) for m in M]); E = np.array(E38)
ZD = np.array([[zeff(dry[m], e) for e in E38] for m in M]); ZW = np.array([[zeff(wet[m], e) for e in E38] for m in M])
fine = [10 ** (1 + 1.0 * i / 399) / 1000.0 for i in range(400)]
def cross(c):
    rf = mass_attenuation(c, fine); d = [p - q for p, q in zip(rf['photoelectric'], rf['incoherent'])]
    for i in range(len(d) - 1):
        if d[i] > 0 >= d[i + 1]: t = d[i] / (d[i] - d[i + 1]); return 1000 * math.exp(math.log(fine[i]) + t * (math.log(fine[i + 1]) - math.log(fine[i])))
CD = np.array([cross(dry[m]) for m in M]); CW = np.array([cross(wet[m]) for m in M])
def ebf(c, e, x):
    z, z1, z2 = zeq(c, e); b, cc, a, Xk, d = material_gp_coeffs(e, z, z1, z2, 'exposure'); return gp_buildup_factor(x, b, cc, a, Xk, d), z
E37 = [e for e in E38 if e >= 0.015]
BD = {e: [ebf(dry[m], e, 40) for m in M] for e in E37}; BW = {e: [ebf(wet[m], e, 40) for m in M] for e in E37}
ix = lambda e: int(np.argmin(abs(E - e)))
sp = lambda v: 100 * (v.max() / v.min() - 1)
out = []; P = lambda *a: print(*a)
P('== MAC change wet vs dry (%), range over the ten mixtures')
for e in (0.01, 0.03, 0.1, 1, 10, 15):
    j = ix(e); d = 100 * (MW[:, j] / MD[:, j] - 1); P('  E=%-5g %+.2f to %+.2f (water MAC %.4g)' % (e, d.min(), d.max(), mass_attenuation(WATER, [e])['total_with_coherent'][0]))
P('== inter-mixture MAC spread (XCOM) dry -> wet, and (max,min) mixtures')
for e in (0.01, 0.03, 0.1, 1, 10, 15):
    j = ix(e); P('  E=%-5g dry %.2f%% (max Mix%d min Mix%d) -> wet %.2f%% (max Mix%d min Mix%d)' % (e, sp(MD[:, j]), MD[:, j].argmax() + 1, MD[:, j].argmin() + 1, sp(MW[:, j]), MW[:, j].argmax() + 1, MW[:, j].argmin() + 1))
P('== Zeff range dry -> wet, (max,min)')
for e in (0.01, 0.02, 0.03, 0.1, 1, 15):
    j = ix(e); P('  E=%-5g dry %.3f-%.3f (max Mix%d min Mix%d) -> wet %.3f-%.3f (max Mix%d min Mix%d); water %.3f' % (e, ZD[:, j].min(), ZD[:, j].max(), ZD[:, j].argmax() + 1, ZD[:, j].argmin() + 1, ZW[:, j].min(), ZW[:, j].max(), ZW[:, j].argmax() + 1, ZW[:, j].argmin() + 1, zeff(WATER, e)))
P('== crossover keV dry %.2f-%.2f -> wet %.2f-%.2f (Mix1 %.2f -> %.2f)' % (CD.min(), CD.max(), CW.min(), CW.max(), CD[0], CW[0]))
P('== density: dry %.3f-%.3f -> wet(additive) %.3f-%.3f ; order same: %s ; highest Mix%d lowest Mix%d' % (min(rho.values()), max(rho.values()), min(rho_wet.values()), max(rho_wet.values()), list(np.argsort([rho[m] for m in M])) == list(np.argsort([rho_wet[m] for m in M])), max(M, key=lambda m: rho_wet[m]), min(M, key=lambda m: rho_wet[m])))
RD = np.array([rho[m] for m in M]); RW = np.array([rho_wet[m] for m in M])
lsp = lambda mac_, r_: np.array([sp(mac_[:, j] * r_) for j in range(len(E))])
a, b_, c_ = lsp(MD, RD), lsp(MW, RW), lsp(MW, RD)
P('== LAC spread over the grid: dry %.1f-%.1f%% | wet additive density %.1f-%.1f%% | wet unchanged density %.1f-%.1f%%' % (a.min(), a.max(), b_.min(), b_.max(), c_.min(), c_.max()))
for e in (0.01, 0.1, 1, 15): j = ix(e); P('     E=%-5g LAC spread dry %.1f | wet additive %.1f | wet same density %.1f' % (e, a[j], b_[j], c_[j]))
j = ix(0.1)
for lab, m_, r_ in (('dry', MD, RD), ('wet additive', MW, RW), ('wet same density', MW, RD)):
    lac = m_[:, j] * r_; P('  at 0.1 MeV %-17s Mix10 needs %.1f%% less thickness than Mix2 ; areal mass Mix10 vs Mix2 %+.2f%% ; largest LAC Mix%d smallest Mix%d' % (lab, 100 * (1 - lac[1] / lac[9]), 100 * (m_[1, j] / m_[9, j] - 1), lac.argmax() + 1, lac.argmin() + 1))
P('== LAC ordering: largest LAC mixture count (of 38) dry/wet-additive:', int(sum((MD[:, k] * RD).argmax() == 9 for k in range(len(E)))), int(sum((MW[:, k] * RW).argmax() == 9 for k in range(len(E)))), '; smallest==Mix2:', int(sum((MD[:, k] * RD).argmin() == 1 for k in range(len(E)))), int(sum((MW[:, k] * RW).argmin() == 1 for k in range(len(E)))))
P('== EBF at 40 mfp, dry -> wet, and Zeq')
for e in (0.04, 0.1, 1, 15):
    d = np.array([v[0] for v in BD[e]]); w = np.array([v[0] for v in BW[e]]); zd = np.array([v[1] for v in BD[e]]); zw = np.array([v[1] for v in BW[e]])
    P('  E=%-5g spread dry %.1f%% (hi Mix%d lo Mix%d) -> wet %.1f%% (hi Mix%d lo Mix%d); EBF change %+.1f%% to %+.1f%%; Zeq order same as EBF order (wet): %s' % (e, sp(d), d.argmax() + 1, d.argmin() + 1, sp(w), w.argmax() + 1, w.argmin() + 1, (100 * (w / d - 1)).min(), (100 * (w / d - 1)).max(), list(np.argsort(zw)) == list(np.argsort(-w))))
# elemental contrast on the wet mean composition
MU = {el: mass_attenuation({el: 1.0}, E38)['total_with_coherent'] for el in 'CHNOS'}
wm_d = {el: np.mean([dry[m].get(el, 0) for m in M]) for el in 'CHNOS'}; wm_w = {el: np.mean([wet[m].get(el, 0) for m in M]) for el in 'CHNOS'}
P('== elemental contrast: % MAC change per +1 wt% replacing C (mean composition), dry -> wet')
for e in (0.01, 0.03, 0.1, 1):
    j = ix(e); bd = sum(wm_d[el] * MU[el][j] for el in 'CHNOS'); bw = sum(wm_w[el] * MU[el][j] for el in 'CHNOS')
    P('  E=%-5g ' % e + ' '.join('%s %+.2f->%+.2f' % (el, 100 * 0.01 * (MU[el][j] - MU['C'][j]) / bd, 100 * 0.01 * (MU[el][j] - MU['C'][j]) / bw) for el in 'HNOS'))
with open(O + 'moisture_mac_zeff.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['E_MeV'] + [f'MACdry_Mix{m}' for m in M] + [f'MACwet_Mix{m}' for m in M] + [f'Zdry_Mix{m}' for m in M] + [f'Zwet_Mix{m}' for m in M])
    for j, e in enumerate(E38): w.writerow([e] + ['%.6g' % x for x in MD[:, j]] + ['%.6g' % x for x in MW[:, j]] + ['%.5g' % x for x in ZD[:, j]] + ['%.5g' % x for x in ZW[:, j]])
with open(O + 'moisture_ebf40.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['E_MeV'] + [f'EBF40dry_Mix{m}' for m in M] + [f'EBF40wet_Mix{m}' for m in M])
    for e in E37: w.writerow([e] + ['%.6g' % v[0] for v in BD[e]] + ['%.6g' % v[0] for v in BW[e]])
json.dump({'water_fraction': FW, 'water_composition': WATER, 'rho_water': RHOW, 'rho_wet_additive': rho_wet}, open(O + 'moisture_assumptions.json', 'w'), indent=1)

# ------------------------------------------------------------------ broad-beam comparison, dry vs wet (same method as analysis/broadbeam/equal_thickness.py)
from scipy import stats
def broad(kind, val, wetflag):
    comps = wet if wetflag else dry; macs = MW if wetflag else MD; rr = RW if wetflag else RD; res = []
    for e in E37:
        j = ix(e); d = macs[:, j] * rr * val if kind == 'thickness' else macs[:, j] * val
        if not np.all((d >= 0.5) & (d <= 40)): continue
        Tp = np.exp(-d); Tb = np.array([ebf(comps[m], e, d[k])[0] for k, m in enumerate(M)]) * Tp
        res.append((e, sp(Tp), sp(Tb), int(Tp.argmin()) + 1, int(Tb.argmin()) + 1, stats.spearmanr(Tp, Tb)[0]))
    return res
TOK = {}
for kind, val, tag in (('thickness', 10.0, 'T'), ('mass', 20.0, 'M')):
    for wetflag, wt in ((False, 'd'), (True, 'w')):
        r = broad(kind, val, wetflag); dd = {round(x[0], 6): x for x in r}
        TOK[f'bb{tag}{wt}_n'] = len(r); TOK[f'bb{tag}{wt}_diff'] = sum(x[3] != x[4] for x in r)
        for e, k in ((0.1, '01'), (1.0, '1')):
            if round(e, 6) in dd: TOK[f'bb{tag}{wt}_p{k}'] = round(dd[round(e, 6)][1], 1); TOK[f'bb{tag}{wt}_t{k}'] = round(dd[round(e, 6)][2], 1)
        band = [x for x in r if 0.1 <= x[0] <= 1.0]; rat = [x[1] / x[2] for x in band]; TOK[f'bb{tag}{wt}_rlo'] = round(min(rat), 1); TOK[f'bb{tag}{wt}_rhi'] = round(max(rat), 1)
        TOK[f'bb{tag}{wt}_spmed'] = round(float(np.median([x[5] for x in r])), 2); TOK[f'bb{tag}{wt}_Emax'] = r[-1][0]
        if round(0.1, 6) in dd: TOK[f'bb{tag}{wt}_bestp'] = dd[0.1][3]; TOK[f'bb{tag}{wt}_bestb'] = dd[0.1][4]
P('== broad-beam dry vs wet'); [P('  ', k, v) for k, v in TOK.items()]
# tokens for the text
def rng2(a, f='%.1f'): return f % a.min() + '--' + f % a.max()
for e, k in ((0.01, '001'), (0.03, '003'), (0.1, '010'), (1, '100'), (10, '1000'), (15, '1500')):
    j = ix(e); TOK['dmac' + k] = rng2(100 * (MW[:, j] / MD[:, j] - 1)); TOK['spd' + k] = '%.2f' % sp(MD[:, j]); TOK['spw' + k] = '%.2f' % sp(MW[:, j])
for e, k in ((0.01, '001'), (1, '100'), (15, '1500')):
    j = ix(e); TOK['zd' + k] = '%.2f--%.2f' % (ZD[:, j].min(), ZD[:, j].max()); TOK['zw' + k] = '%.2f--%.2f' % (ZW[:, j].min(), ZW[:, j].max())
TOK['zmax4_wet'] = bool(all(ZW[:, j].argmax() == 3 for j in range(len(E)))); TOK['zmin2_wet'] = bool(all(ZW[:, j].argmin() == 1 for j in range(len(E)) if E[j] >= 0.02)); TOK['zmin10_wet_001'] = bool(ZW[:, 0].argmin() == 9)
TOK['crd'] = '%.1f--%.1f' % (CD.min(), CD.max()); TOK['crw'] = '%.1f--%.1f' % (CW.min(), CW.max())
TOK['rhod'] = '%.3f--%.3f' % (RD.min(), RD.max()); TOK['rhow'] = '%.3f--%.3f' % (RW.min(), RW.max()); TOK['rhoorder'] = bool(list(np.argsort(RD)) == list(np.argsort(RW)))
TOK['lspd'] = '%.1f--%.1f' % (a.min(), a.max()); TOK['lspw'] = '%.1f--%.1f' % (b_.min(), b_.max()); TOK['lspw2'] = '%.1f--%.1f' % (c_.min(), c_.max())
for lab, m_, r_ in (('d', MD, RD), ('w', MW, RW), ('s', MW, RD)):
    lac = m_[:, j0] * r_ if False else m_[:, ix(0.1)] * r_; TOK['thin' + lab] = '%.1f' % (100 * (1 - lac[1] / lac[9])); TOK['am' + lab] = '%.2f' % (100 * (m_[1, ix(0.1)] / m_[9, ix(0.1)] - 1))
for e, k in ((0.04, '004'), (0.1, '010'), (1, '100')):
    d_ = np.array([v[0] for v in BD[e]]); w_ = np.array([v[0] for v in BW[e]])
    TOK['ebfd' + k] = '%.1f' % sp(d_); TOK['ebfw' + k] = '%.1f' % sp(w_); TOK['ebfdhi' + k] = int(d_.argmax()) + 1; TOK['ebfwhi' + k] = int(w_.argmax()) + 1; TOK['ebfdlo' + k] = int(d_.argmin()) + 1; TOK['ebfwlo' + k] = int(w_.argmin()) + 1
    TOK['ebfchg' + k] = '%.0f--%.0f' % tuple(sorted((-(100 * (w_ / d_ - 1)).max(), -(100 * (w_ / d_ - 1)).min())))
zw10 = np.array([v[1] for v in BW[0.1]]); TOK['zeqmonowet'] = bool(list(np.argsort(zw10)) == list(np.argsort(-np.array([v[0] for v in BW[0.1]]))))
for el in 'HNOS':
    j = ix(0.01); bd = sum(wm_d[x] * MU[x][j] for x in 'CHNOS'); bw = sum(wm_w[x] * MU[x][j] for x in 'CHNOS'); TOK['cd' + el] = '%+.2f' % (100 * 0.01 * (MU[el][j] - MU['C'][j]) / bd); TOK['cw' + el] = '%+.2f' % (100 * 0.01 * (MU[el][j] - MU['C'][j]) / bw)
json.dump(TOK, open(O + 'moisture_tokens.json', 'w'), indent=1); print('tokens', len(TOK))
