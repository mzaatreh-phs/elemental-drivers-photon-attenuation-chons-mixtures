#!/usr/bin/env python3
"""Rewrite Results 3.2-3.7 of CHON_shielding(20).tex.  Every number in the prose is computed here from the analysis outputs (tokens @name@)."""
import os, sys, csv, re, json, math
import numpy as np
A = os.path.expanduser('~/CHON/recalc_2026-09-19/analysis/')
sys.path.insert(0, A + 'g4')
from g4analysis import *
TEX = os.path.expanduser('~/CHON/CHOM_submitted version/CHON_shielding(20).tex')
s = open(TEX).read()
def rows(p): return list(csv.DictReader(open(A + p)))

# ------------------------------------------------------------------ data
runs, _ = load_runs(); en = complete_energies(runs); D = build_tables(runs, en); E = np.array(D['E'])
xc = load_xcom_mac(); X = np.array([[xc[(m, round(e, 6))] for e in E] for m in MIXES])
mac, lac, u, T = D['mac'], D['lac'], D['u'], D['T']; rho = np.array([RHO[m] for m in MIXES]); Wd = {el: np.array([COMP[m].get(el, 0.0) for m in MIXES]) for el in 'CHONS'}
ix = lambda e: int(np.argmin(abs(E - e)))
C5 = rows('g4/results/table5_correlations.csv')
def c5(p, r, e, k='r_G4'): return float([x for x in C5 if x['pred'] == p and x['resp'] == r and abs(float(x['E']) - e) < 1e-9][0][k])
T4 = {t['E']: t for t in json.load(open(A + 'g4/results/table4_spreads.json'))}
V = {}
def sg(x): return '%+.2f' % x
def put(k, v, f='%.2f'): V[k] = f % v if not isinstance(v, str) else v
# --- library coupling
for k, (a, b) in {'rOC': ('O', 'C'), 'rON': ('O', 'N'), 'rOS': ('O', 'S')}.items(): V[k] = sg(np.corrcoef(Wd[a], Wd[b])[0, 1])
V['rOD'] = sg(np.corrcoef(Wd['O'], rho)[0, 1])
# --- Table 5 values used in the text
for tag, e in (('001', 0.01), ('003', 0.03), ('005', 0.05), ('010', 0.1), ('100', 1.0), ('1000', 10.0), ('1500', 15.0)):
    V['rN' + tag] = sg(c5('Nitrogen', 'MAC', e)); V['rO' + tag] = sg(c5('Oxygen', 'MAC', e)); V['rH' + tag] = sg(c5('Hydrogen', 'MAC', e))
    V['xO' + tag] = sg(c5('Oxygen', 'MAC', e, 'r_xcom')); V['xH' + tag] = sg(c5('Hydrogen', 'MAC', e, 'r_xcom'))
    V['rDM' + tag] = sg(c5('Density', 'MAC', e)); V['rDL' + tag] = '%+.3f' % c5('Density', 'LAC', e)
    V['rSZ' + tag] = sg(c5('Sulfur', 'Zeff', e)); V['rOZ' + tag] = sg(c5('Oxygen', 'Zeff', e))
# --- spreads
xsp = lambda e: spread_pct(X[:, ix(e)]); gsp = lambda e: spread_pct(mac[:, ix(e)])
for tag, e in (('001', 0.01), ('003', 0.03), ('010', 0.1), ('100', 1), ('10', 10), ('15', 15)): put('sp' + tag, gsp(e)); put('xsp' + tag, xsp(e))
V['mac001lo'] = '%.2f' % mac[:, 0].min(); V['mac001hi'] = '%.2f' % mac[:, 0].max()
V['band10'] = '%.2f--%.2f' % (T4[10.0]['MAC']['mc_16'], T4[10.0]['MAC']['mc_84'])
# --- process fractions / crossover
P = [r for r in rows('xcom/process_fractions.csv') if r['set'] == 'new' and r['material'] == 'Mix1']; Pd = {round(float(r['E_MeV']), 6): r for r in P}
V['pe001'] = '%.1f' % (100 * float(Pd[0.01]['F_photoelectric'])); V['pe010'] = '%.1f' % (100 * float(Pd[0.1]['F_photoelectric']))
V['inc001'] = '%.1f' % (100 * float(Pd[0.01]['F_incoherent'])); V['inc010'] = '%.1f' % (100 * float(Pd[0.1]['F_incoherent'])); V['pair15'] = '%.0f' % (100 * float(Pd[15.0]['F_pair']))
V['coh_pk'] = '%.1f' % max(100 * float(r['F_coherent']) for r in P)
V['coh_max'] = '%.1f' % max(100 * float(r['F_coherent']) for r in rows('xcom/process_fractions.csv') if r['set'] == 'new')
CR = {r['material']: float(r['crossover_keV']) for r in rows('xcom/crossover.csv') if r['set'] == 'new'}
V['cr_lo'] = '%.1f' % min(CR.values()); V['cr_hi'] = '%.1f' % max(CR.values()); V['cr1'] = '%.1f' % CR['Mix1']
# --- elemental sensitivity table + text values
SENS = {float(r['Energy_MeV']): r for r in rows('sensitivity/exchange_sensitivity_and_contribution.csv')}
ENS = {float(r['Energy_MeV']): r for r in rows('sensitivity/ensemble_vs_library_correlation_and_importance.csv')}
def dm(e, el): return float(SENS[e]['pct_MAC_change_per_1wt%%_%s_replacing_C' % el])
for tag, e in (('001', 0.01), ('003', 0.03), ('005', 0.05), ('010', 0.1), ('100', 1.0)):
    for el in 'HNOS': V['d%s%s' % (el, tag)] = '%+.2f' % dm(e, el) if abs(dm(e, el)) >= 0.01 else '%+.3f' % dm(e, el)
V['dS001p'] = '%.1f' % dm(0.01, 'S'); V['dS003p'] = '%.1f' % dm(0.03, 'S'); V['dS005p'] = '%.1f' % dm(0.05, 'S'); V['dS010p'] = '%.2f' % dm(0.1, 'S'); V['dS100p'] = '%.3f' % dm(1.0, 'S')
V['dO001p'] = '%.2f' % dm(0.01, 'O'); V['dN001p'] = '%.2f' % dm(0.01, 'N'); V['dH001p'] = '%.2f' % dm(0.01, 'H')
V['ratioSO'] = '%.0f' % (dm(0.01, 'S') / dm(0.01, 'O')); V['ratioSN'] = '%.0f' % (dm(0.01, 'S') / dm(0.01, 'N'))
sel = [e for e in SENS if 0.1 <= e <= 3.0]; V['dHmin'] = '%.2f' % min(dm(e, 'H') for e in sel); V['dHmax'] = '%.2f' % max(dm(e, 'H') for e in sel)
V['dONmax'] = '%.2f' % (math.ceil(100 * max(max(abs(dm(e, 'O')), abs(dm(e, 'N'))) for e in sel)) / 100)
def imp(e, el): return float(ENS[e]['importance_pct_' + el])
V['ensO001'] = sg(float(ENS[0.01]['r_ensemble_O'])); V['libO001'] = sg(float(ENS[0.01]['r_library_O'])); V['ensS001'] = sg(float(ENS[0.01]['r_ensemble_S'])); V['libS001'] = sg(float(ENS[0.01]['r_library_S']))
for el in 'HNOS': V['imp%s001' % el] = '%.0f' % imp(0.01, el)
V['impH100'] = '%.0f' % imp(1.0, 'H')
RT = rows('sensitivity/range_test_importance.csv'); scen = sorted({r['scenario'] for r in RT})
def dom(e):
    out = {}
    for sc in scen:
        r = [x for x in RT if x['scenario'] == sc and abs(float(x['Energy_MeV']) - e) < 1e-9]
        if r: out[sc] = max('HNOS', key=lambda el: float(r[0]['share_' + el]))
    return out
V['nscen'] = str(len(scen)); en_rt = sorted({float(r['Energy_MeV']) for r in RT})
V['Hall'] = 'all' if all(sum(v == 'H' for v in dom(e).values()) == len(scen) for e in en_rt if 0.1 <= e <= 3.0) else 'NOT-ALL'
V['S001'] = str(sum(v == 'S' for v in dom(0.01).values())); V['S003'] = str(sum(v == 'S' for v in dom(0.03).values())); V['S005'] = str(sum(v == 'S' for v in dom(0.05).values()))
V['H10'] = str(sum(v == 'H' for v in dom(10.0).values())); V['O10'] = str(sum(v == 'O' for v in dom(10.0).values())); V['O15'] = str(sum(v == 'O' for v in dom(15.0).values())); V['H15'] = str(sum(v == 'H' for v in dom(15.0).values()))
EX = [r for r in rows('sensitivity/exchange_confirmation.csv') if r['variant'] in ('dO', 'dN', 'dS', 'dH', 'dALL')]; z = np.array([float(r['z']) for r in EX])
V['zn'] = str(len(z)); V['zmean'] = '%.2f' % z.mean(); V['zstd'] = '%.2f' % z.std(ddof=1); V['zchi'] = '%.2f' % (z ** 2).mean(); V['zmax'] = '%.1f' % abs(z).max()
def ex(e, v): r = [x for x in EX if abs(float(x['E_MeV']) - e) < 1e-9 and x['variant'] == v][0]; return float(r['dMAC_Geant4_pct']), float(r['u_pct']), float(r['dMAC_XCOM_mixture_rule_pct'])
g, uu, xx = ex(0.01, 'dS'); V['exS_g'] = '%+.1f' % g; V['exS_u'] = '%.1f' % uu; V['exS_x'] = '%+.1f' % xx
g, uu, xx = ex(0.01, 'dO'); V['exO_g'] = '%+.1f' % g; V['exO_u'] = '%.1f' % uu; V['exO_x'] = '%+.1f' % xx
g, uu, xx = ex(0.01, 'dALL'); V['exA_g'] = '%+.1f' % g; V['exA_u'] = '%.1f' % uu; V['exA_x'] = '%+.1f' % xx
# --- LAC / density
V['lacsp001'] = '%.1f' % spread_pct(lac[:, 0]); V['lacsp003'] = '%.1f' % spread_pct(lac[:, ix(0.03)]); V['lacsp010'] = '%.1f' % spread_pct(lac[:, ix(0.1)])
V['lacsp100'] = '%.1f' % spread_pct(lac[:, ix(1)]); V['lacsp1000'] = '%.1f' % spread_pct(lac[:, ix(10)]); V['lacsp1500'] = '%.1f' % spread_pct(lac[:, ix(15)])
lsp = np.array([spread_pct(lac[:, j]) for j in range(len(E))]); V['lacspmin'] = '%.1f' % lsp.min(); V['lacspmax'] = '%.1f' % lsp.max()
rl = np.array([pearson(list(rho), list(lac[:, j])) for j in range(len(E))]); V['rLmin'] = '%.3f' % rl.min(); V['rLmax'] = '%.3f' % rl.max()
rm = np.array([pearson(list(rho), list(1 / lac[:, j])) for j in range(len(E))]); V['rMmin'] = '%.3f' % abs(rm).min(); V['rMmax'] = '%.3f' % abs(rm).max()
rr = np.array([pearson(list(rho), list(100 * (1 - T[:, j]))) for j in range(len(E))]); V['rRmin'] = '%.3f' % rr.min(); V['rRmax'] = '%.3f' % rr.max()
nE = len(E); V['nE'] = str(nE); V['n10'] = str(int(sum(lac[:, j].argmax() == 9 for j in range(nE)))); V['n2'] = str(int(sum(lac[:, j].argmin() == 1 for j in range(nE))))
V['n10rpe'] = str(int(sum(T[:, j].argmin() == 9 for j in range(nE)))); V['n2rpe'] = str(int(sum(T[:, j].argmax() == 1 for j in range(nE))))
not10 = [float(E[j]) for j in range(nE) if lac[:, j].argmax() != 9]; V['not10'] = ', '.join('%g' % e for e in not10)
# --- MAC / LAC at selected energies
def rg(a, f): return f % a.min() + '--' + f % a.max()
for tag, e in (('001', 0.01), ('003', 0.03), ('010', 0.1), ('100', 1), ('1000', 10), ('1500', 15)):
    j = ix(e); V['MAC' + tag] = rg(mac[:, j], '%.4g' if e >= 0.1 else '%.3g'); V['LAC' + tag] = rg(lac[:, j], '%.4g'); V['HVL' + tag] = rg(np.log(2) / lac[:, j], '%.3g'); V['MFP' + tag] = rg(1 / lac[:, j], '%.3g'); V['TVL' + tag] = rg(np.log(10) / lac[:, j], '%.3g')
    V['RPE' + tag] = rg(100 * (1 - T[:, j]), '%.2f')
V['MAC001'] = '%.2f--%.2f' % (mac[:, 0].min(), mac[:, 0].max()); V['LAC001'] = '%.2f--%.2f' % (lac[:, 0].min(), lac[:, 0].max())
for tag, e in (('010', 0.1), ('1000', 10)):
    j = ix(e); V['MAC' + tag] = '%.4f--%.4f' % (mac[:, j].min(), mac[:, j].max()) if e < 1 else '%.5f--%.5f' % (mac[:, j].min(), mac[:, j].max())
    V['LAC' + tag] = '%.4f--%.4f' % (lac[:, j].min(), lac[:, j].max()) if e < 1 else '%.5f--%.5f' % (lac[:, j].min(), lac[:, j].max())
j = ix(0.1); l10, l2, m10, m2 = lac[9, j], lac[1, j], mac[9, j], mac[1, j]; um10, um2 = u[9, j], u[1, j]
V['thin'] = '%.1f' % (100 * (1 - l2 / l10)); V['am'] = '%.1f' % (100 * (m2 / m10 - 1)); V['amu'] = '%.1f' % (100 * (m2 / m10) * math.hypot(um2 / m2, um10 / m10)); V['amx'] = '%.1f' % (100 * (X[1, j] / X[9, j] - 1))
V['H2'] = '%.4f' % Wd['H'][1]; V['H10'] = V.get('H10', ''); V['Hm10'] = '%.4f' % Wd['H'][9]
# --- Zeff
Z = {}
for r in rows('xcom/zeff.csv'):
    if r['set'] in ('new', 'ref'): Z[(r['material'], round(float(r['E_MeV']), 6))] = float(r['Zeff_direct'])
def zb(e): v = np.array([Z[('Mix%d' % m, e)] for m in range(1, 11)]); return v
for tag, e in (('001', 0.01), ('020', 0.02), ('030', 0.03), ('100', 1.0), ('1000', 10.0), ('1500', 15.0)):
    v = zb(e); V['Zlo' + tag] = '%.2f' % v.min(); V['Zhi' + tag] = '%.2f' % v.max(); V['Zsp' + tag] = '%.2f' % (v.max() - v.min()); V['Zwat' + tag] = '%.2f' % Z[('water', e)]; V['Zpm' + tag] = '%.2f' % Z[('PMMA', e)]
V['Z1lo'] = '%.2f' % zb(1.0).min(); V['Z1hi'] = '%.2f' % zb(1.0).max()
Zmax4 = all(zb(e).argmax() == 3 for e in sorted({k[1] for k in Z if k[0] == 'Mix1'})); V['Zmax4'] = 'yes' if Zmax4 else 'NO'
Zmin2 = all(zb(e).argmin() == 1 for e in sorted({k[1] for k in Z if k[0] == 'Mix1'}) if e >= 0.02); V['Zmin2'] = 'yes' if Zmin2 else 'NO'
V['Zwat_first_below'] = '%g' % min(e for e in sorted({k[1] for k in Z if k[0] == 'water'}) if Z[('water', e)] < zb(e).min() - 1e-9)
V['Hwater'] = '11.2'; V['Hmean'] = '%.1f' % (100 * Wd['H'].mean()); V['Hlo'] = '%.1f' % (100 * Wd['H'].min()); V['Hhi'] = '%.1f' % (100 * Wd['H'].max())
V['z10_4'] = '%.3f' % Z[('Mix4', 10.0)]; V['z10_10'] = '%.3f' % Z[('Mix10', 10.0)]
# --- EBF
EB = [r for r in rows('xcom/exposure_buildup.csv') if r['set'] == 'new']; Es = sorted({float(r['E_MeV']) for r in EB})
def Bf(m, e, d): return float([r for r in EB if r['material'] == 'Mix%d' % m and abs(float(r['E_MeV']) - e) < 1e-9][0]['B_%dmfp' % d])
def Bv(e, d): return np.array([Bf(m, e, d) for m in range(1, 11)])
def pk(d): return [max(Bf(m, e, d) for e in Es) for m in range(1, 11)]
V['pk5'] = '%.0f--%.0f' % (min(pk(5)), max(pk(5))); V['pk10'] = '%.0f--%.0f' % (min(pk(10)), max(pk(10))); V['pk40lo'] = '%.1f' % (min(pk(40)) / 1e4); V['pk40hi'] = '%.1f' % (max(pk(40)) / 1e4)
for tag, e in (('0015', 0.015), ('002', 0.02), ('003', 0.03)): v = Bv(e, 40); V['E40' + tag] = '%.2f--%.2f' % (v.min(), v.max()) if v.max() < 10 else '%.1f--%.1f' % (v.min(), v.max())
v4 = Bv(0.04, 40); V['ebfsp004'] = '%.1f' % (100 * (v4.max() / v4.min() - 1)); V['b004lo'] = '%.0f' % v4.min(); V['b004hi'] = '%.0f' % v4.max(); V['m004lo'] = str(v4.argmin() + 1); V['m004hi'] = str(v4.argmax() + 1)
v1 = Bv(0.1, 40); V['ebfsp010'] = '%.1f' % (100 * (v1.max() / v1.min() - 1)); V['b010lo'] = '{:,.0f}'.format(v1.min()).replace(',', '{,}'); V['b010hi'] = '{:,.0f}'.format(v1.max()).replace(',', '{,}'); V['m010lo'] = str(v1.argmin() + 1); V['m010hi'] = str(v1.argmax() + 1)
V['red010'] = '%.1f' % (100 * (1 - v1.min() / v1.max()))
V['proxlo'] = '%.2f' % (1e13 * math.exp(-40) * v1.min()); V['proxhi'] = '%.2f' % (1e13 * math.exp(-40) * v1.max())
v15 = Bv(15.0, 40); V['b15lo'] = '%.3f' % v15.min(); V['b15hi'] = '%.3f' % v15.max(); V['m15lo'] = str(v15.argmin() + 1); V['m15hi'] = str(v15.argmax() + 1); V['sp15e'] = '%.1f' % (100 * (v15.max() / v15.min() - 1))
v1m = Bv(1.0, 40); V['sp1e'] = '%.1f' % (100 * (v1m.max() / v1m.min() - 1)); v10 = Bv(10.0, 40); V['sp10e'] = '%.1f' % (100 * (v10.max() / v10.min() - 1))
Zq = {int(r['material'][3:]): float(r['Zeq']) for r in EB if abs(float(r['E_MeV']) - 0.1) < 1e-9}; V['Zq10'] = '%.2f' % Zq[10]; V['Zq4'] = '%.2f' % Zq[4]
mono = [m for m, _ in sorted(Zq.items(), key=lambda t: t[1])]; ebf_order = [m for m in sorted(range(1, 11), key=lambda m: -Bf(m, 0.1, 40))]; V['zeqmono'] = 'yes' if mono == ebf_order else 'NO'; assert mono == ebf_order, 'EBF not monotonic in Zeq'
# --- broad-beam
BS = rows('broadbeam/broadbeam_summary.csv')
def bb(kind, val): return [r for r in BS if r['condition'] == kind and r['value'] == val]
for tag, (kind, val) in {'T': ('thickness_cm', '10'), 'M': ('areal_mass_g_cm2', '20')}.items():
    sub = bb(kind, val); V['bbn' + tag] = str(len(sub)); V['bbdiff' + tag] = str(sum(r['best_mix_primary'] != r['best_mix_broad'] for r in sub))
    d = {round(float(r['E_MeV']), 6): r for r in sub}
    for tg, e in (('015', 0.15), ('01', 0.1), ('1', 1.0)):
        if e in d: V['bbp%s%s' % (tag, tg)] = '%.1f' % float(d[e]['spread_primary_pct']); V['bbt%s%s' % (tag, tg)] = '%.1f' % float(d[e]['spread_broad_pct'])
    band = [r for r in sub if 0.1 <= float(r['E_MeV']) <= 1.0]; rat = [float(r['spread_primary_pct']) / float(r['spread_broad_pct']) for r in band]; V['bbr%slo' % tag] = '%.1f' % min(rat); V['bbr%shi' % tag] = '%.1f' % max(rat)
    sp = [float(r['spearman_primary_vs_broad']) for r in sub]; V['bbsp%smed' % tag] = '%.2f' % np.median(sp); V['bbsp%smin' % tag] = '%.2f' % min(sp)
    hi = [r for r in sub if float(r['E_MeV']) >= 1.3]; V['bbhi%s' % tag] = '%.2f' % max(float(r['spread_broad_pct']) for r in hi); V['bbhip%s' % tag] = '%.2f' % max(float(r['spread_primary_pct']) for r in hi)
    V['bbEmax' + tag] = '%g' % float(sub[-1]['E_MeV']); V['bbEmin' + tag] = '%g' % float(sub[0]['E_MeV'])
dT = {round(float(r['E_MeV']), 6): r for r in bb('thickness_cm', '10')}; V['bbT_bestp'] = dT[0.1]['best_mix_primary']; V['bbT_bestb'] = dT[0.1]['best_mix_broad']
dM = {round(float(r['E_MeV']), 6): r for r in bb('areal_mass_g_cm2', '20')}; V['bbM_bestp'] = dM[0.1]['best_mix_primary']; V['bbM_bestb'] = dM[0.1]['best_mix_broad']
V['bbT_1'] = '%g' % min(float(r['E_MeV']) for r in bb('thickness_cm', '10') if r['best_mix_primary'] != r['best_mix_broad'])
# ------------------------------------------------------------------ Table 7 (sensitivity + ensemble importance)
TE = (0.01, 0.03, 0.05, 0.1, 1.0, 10.0, 15.0)
def cell(e, el): v = dm(e, el); return '%+.2f' % v if abs(v) >= 0.005 else '%+.3f' % v
tab7 = (r'''\begin{table}[htbp]
\centering
\caption{Elemental sensitivity of the MAC at selected photon energies (XCOM mixture rule, base composition equal to the mean of the ten mixtures): change in MAC (\%) per +1~wt\% of the element replacing carbon (left), and importance (\%) of each element for the MAC variability of the decorrelated ensemble (right; Section~\ref{sec:sensitivity}). The importance ranking depends on the assumed element ranges (see text).}
\label{tbl:sensitivity}
\small
\setlength{\tabcolsep}{5pt}
\begin{tabular}{c cccc cccc}
\toprule
 & \multicolumn{4}{c}{$\Delta\mu_m/\mu_m$ per +1~wt\% (\%)} & \multicolumn{4}{c}{Importance (\%)} \\
\cmidrule(lr){2-5}\cmidrule(lr){6-9}
$E$ (MeV) & H & N & O & S & H & N & O & S \\
\midrule
''' + '\n'.join('%-5g & ' % e + ' & '.join(['$%s$' % cell(e, el) for el in 'HNOS']) + ' & ' + ' & '.join('%.0f' % imp(e, el) for el in 'HNOS') + r' \\' for e in TE) + r'''
\bottomrule
\end{tabular}
\end{table}''')

# ------------------------------------------------------------------ keep existing float environments
def grab(label, kind):
    i = s.index(r'\label{%s}' % label); a = s.rindex(r'\begin{%s}' % kind, 0, i); b = s.index(r'\end{%s}' % kind, i) + len(r'\end{%s}' % kind); return s[a:b]
figP = grab('fig:process', 'figure'); tabM = grab('tbl:macspread', 'table'); tabC = grab('tbl:correlation', 'table')
figML = grab('fig:maclac', 'figure'); figMF = grab('fig:mfphvltvl', 'figure'); figR = grab('fig:rpe', 'figure'); figZ = grab('fig:zeff', 'figure'); figZB = grab('fig:zeff_baseline', 'figure'); tabZ = grab('tbl:zeff', 'table')
figE = grab('fig:ebf', 'figure'); figBB = grab('fig:broadbeam', 'figure')
figP = re.sub(r'\\caption\{.*?\}\n', r'\\caption{Process-fraction decomposition of the total MAC of Mixture~1 (XCOM). The dashed line marks the photoelectric--incoherent crossover at @cr1@~keV.}' + '\n', figP, flags=re.S)
figML = re.sub(r'\\caption\{.*?\}\n', r'\\caption{Mass attenuation coefficient (left) and linear attenuation coefficient (right) versus photon energy for the ten mixtures. Insets resolve the low-energy region.}' + '\n', figML, flags=re.S)
figMF = re.sub(r'\\caption\{.*?\}\n', r'\\caption{Half-value layer (left), tenth-value layer (middle), and mean free path (right) versus photon energy for the ten mixtures.}' + '\n', figMF, flags=re.S)
figR = re.sub(r'\\caption\{.*?\}\n', r'\\caption{Radiation protection efficiency (primary-beam attenuation only) of the ten mixtures at a fixed thickness of 1~cm.}' + '\n', figR, flags=re.S)
figZ = re.sub(r'\\caption\{.*?\}\n', r'\\caption{Effective atomic number versus photon energy for the ten mixtures (XCOM direct method).}' + '\n', figZ, flags=re.S)
figZB = re.sub(r'\\caption\{.*?\}\n', r'\\caption{Effective atomic number of the ten mixtures (shaded band) compared with water and PMMA.}' + '\n', figZB, flags=re.S)
figE = re.sub(r'\\caption\{.*?\}\n', r'\\caption{Exposure buildup factor (G-P) of the ten mixtures at penetration depths of 5--40~mfp.}' + '\n', figE, flags=re.S)

# ------------------------------------------------------------------ text: 3.2
sec32 = (r'''\subsection{Composition Sensitivity and Elemental Drivers Across Photon-Interaction Regimes}\label{sec:composition}

The ten formulations span oxygen mass fractions of 0.2368--0.4586 and densities of 1.268--1.415~g/cm\textsuperscript{3}, but their elemental fractions are strongly coupled (Section~\ref{sec:materials}): oxygen correlates negatively with carbon ($r=@rOC@$), with nitrogen and with sulfur ($r=@rON@$ for both, which are collinear) and positively with density ($r=@rOD@$). Because MAC is normalized by density, density cannot enter it mechanically, so any density--MAC association is an artefact of the composition coupling. Correlations across the library are therefore read as descriptions of coupled composition trends, and the effects of individual elements are separated afterwards by the elemental-sensitivity analysis.

In the photoelectric region the MAC spread is @sp001@\% at 0.01~MeV and @sp003@\% at 0.03~MeV (Table~\ref{tbl:macspread}), and MAC correlates most strongly and positively with the nitrogen, sulfur, and BSA fractions ($r=@rN001@$ at 0.01~MeV, $@rN003@$ at 0.03~MeV, and $@rN005@$ at 0.05~MeV; Table~\ref{tbl:correlation}), whereas the oxygen coefficient is negative ($@rO001@$, $@rO003@$, and $@rO005@$; noise-free XCOM values $@xO001@$, $@xO003@$, and $@xO005@$). An argument that the heaviest CHON element controls the photoelectric response would predict the opposite, and the result follows from the design: the oxygen-richest mixtures are those with the least BSA, and BSA is the only carrier of sulfur and nitrogen. The hydrogen coefficient changes sign between 0.03 and 0.05~MeV ($@rH001@$ at 0.01~MeV, $@rH003@$ at 0.03~MeV, $@rH005@$ at 0.05~MeV) and reaches $@rH010@$ at 0.1~MeV; the noise-free values are $@xH010@$ at 0.1~MeV and $@xH100@$ at 1~MeV, so that in the Compton region the small residual MAC spread is systematically controlled by hydrogen. Above 1~MeV, however, no element sustains a consistent association (noise-free hydrogen coefficient $@xH1000@$ at 10~MeV and $@xH1500@$ at 15~MeV), and the Geant4 rows of Table~\ref{tbl:correlation} are limited by counting noise there.

Density also correlates with MAC at low energy ($r=@rDM001@$ at 0.01~MeV). That association cannot be causal, since MAC is normalized by density; it reflects the coupling of density to the BSA-related elements, and it weakens and reverses as the photoelectric contribution disappears ($@rDM010@$ at 0.1~MeV), whereas the density--LAC coefficient is $@rDL001@$ at 0.01~MeV and $@rDL010@$ at 0.1~MeV.

The process decomposition of Eqs.~\eqref{eq:processmix} and~\eqref{eq:processfrac}, shown in Figure~\ref{fig:process}, supplies the physical explanation for the energy dependence. The photoelectric and incoherent contributions cross at @cr_lo@--@cr_hi@~keV across the ten mixtures, at @cr1@~keV for Mixture~1. For that mixture, the photoelectric contribution falls from @pe001@\% of the total MAC at 0.01~MeV to @pe010@\% at 0.1~MeV, while the incoherent contribution rises from @inc001@\% to @inc010@\% and remains dominant through the intermediate-energy region. The coherent term is a minor correction that peaks at @coh_pk@\% near 0.03~MeV (up to @coh_max@\% across the mixtures). Pair production becomes non-negligible above a few MeV and contributes @pair15@\% of the total MAC of Mixture~1 at 15~MeV. The process crossover near 26~keV marks the change in dominant interaction mechanism.

''' + figP + r'''

The inter-mixture spread defined by Eq.~\eqref{eq:spread} makes the energy dependence of compositional sensitivity explicit (Table~\ref{tbl:macspread}). At 0.01~MeV the MAC ranges from @mac001lo@ to @mac001hi@~cm\textsuperscript{2}/g, a spread of @sp001@\%, reproduced by the noise-free XCOM value of @xsp001@\%. It falls to @sp003@\% at 0.03~MeV and to about 1\% above 0.08~MeV. From 0.1~MeV upward the identities of the maximum- and minimum-MAC mixtures are not statistically established, and at 10~MeV the Geant4 spread (@sp10@\%) is comparable to its own noise band (@band10@\%) and well above the noise-free XCOM value (@xsp10@\%): the composition sensitivity of mass-normalized attenuation above about 0.1~MeV is small (about 1\%) and cannot be resolved between individual mixtures at $2\times10^{6}$ histories.

''' + tabM + '\n\n' + tabC + r'''

The elemental-sensitivity analysis separates these effects (Table~\ref{tbl:sensitivity}). At 0.01~MeV, replacing 1~wt\% of carbon by sulfur raises the MAC by @dS001p@\%, about @ratioSO@ times the effect of the same mass of oxygen (@dO001p@\%) and @ratioSN@ times that of nitrogen (@dN001p@\%), whereas hydrogen lowers it ($@dH001p@$\%). The sulfur effect falls to @dS003p@\% at 0.03~MeV, @dS005p@\% at 0.05~MeV, and @dS010p@\% at 0.1~MeV, and is negligible near 1~MeV (@dS100p@\%). Hydrogen is the only element whose contrast with carbon changes sign: it lowers the MAC below, and raises it above, about 23--24~keV (the exact energy depends on the interpolation between XCOM grid points), and between 0.1 and 3~MeV it raises the MAC by @dHmin@--@dHmax@\% per weight percent, whereas oxygen or nitrogen change it by at most @dONmax@\% per weight percent over the same range.

''' + tab7 + r'''

The library correlations do not reveal this ordering, because sulfur enters only with BSA and BSA covaries with oxygen. In the decorrelated ensemble the oxygen--MAC coefficient at 0.01~MeV is $@ensO001@$ instead of the library value $@libO001@$, and the sulfur coefficient is $@ensS001@$ instead of $@libS001@$; the importance of sulfur for MAC variability is @impS001@\%, against @impO001@\% for oxygen, @impN001@\% for nitrogen, and @impH001@\% for hydrogen (Table~\ref{tbl:sensitivity}). At 1~MeV hydrogen accounts for @impH100@\% of the variability. This ranking depends on the assumed ranges, and only part of it persists across the seven scenarios of Section~\ref{sec:sensitivity}: hydrogen is the dominant element at every energy between 0.1 and 3~MeV in @Hall@ seven scenarios; below 0.06~MeV sulfur dominates in @S001@, @S003@, and @S005@ of the seven scenarios at 0.01, 0.03, and 0.05~MeV, respectively, while oxygen (or, at 0.05~MeV, hydrogen) takes over in the remaining scenarios, in which the ranges are halved or doubled or the compositions are Dirichlet distributed; and above about 8~MeV no element dominates consistently (hydrogen in @H10@ and oxygen in @O10@ scenarios at 10~MeV; oxygen in @O15@ and hydrogen in @H15@ at 15~MeV). The robust statements are therefore that sulfur has by far the largest contrast per unit mass at low energy, that hydrogen controls the small MAC variation in the Compton region, and that oxygen is not a general control variable.

The mixture-rule predictions were confirmed with Geant4. Over the @zn@ comparisons (five composition changes at eight energies) the normalized differences between the Geant4 and mixture-rule changes have a mean of $@zmean@$, a standard deviation of @zstd@, and $\chi^{2}/N=@zchi@$, with none beyond $@zmax@\sigma$. For example, at 0.01~MeV replacing 1.5~wt\% of carbon by sulfur changes the Geant4 MAC by $@exS_g@\pm@exS_u@\%$ against $@exS_x@$\% predicted, replacing 5~wt\% by oxygen by $@exO_g@\pm@exO_u@\%$ against $@exO_x@$\%, and all four exchanges together by $@exA_g@\pm@exA_u@\%$ against $@exA_x@$\%. The changes are resolvable in Geant4 only below about 0.3~MeV, where they exceed the counting uncertainty of 0.14--0.6\%.

Table~\ref{tbl:correlation} collects the correlation coefficients of density with the derived quantities. The behaviour of LAC is different from that of MAC because $\mu=\rho\mu_m$ carries the density explicitly. Density--LAC correlations range from @rLmin@ to @rLmax@ over the 38-energy grid, and the LAC spread ranges from @lacspmin@ to @lacspmax@\% (@lacsp001@\% at 0.01~MeV, @lacsp010@\% at 0.1~MeV, and @lacsp1500@\% at 15~MeV). Mixture~2, the least dense, has the smallest LAC at all @nE@ energies; the largest LAC belongs to Mixture~10, the densest, at @n10@ of them, and to Mixture~3 at the other eight (@not10@~MeV). The same distinction propagates algebraically to MFP, HVL, TVL, and RPE: density correlates with MFP, HVL, and TVL with absolute coefficients of @rMmin@--@rMmax@ and with RPE between @rRmin@ and @rRmax@ over the grid. These derived-parameter correlations are transformations of the same LAC rather than independent evidence.

''')
json.dump(V, open(A + 'results_text/tokens.json', 'w'), indent=1)

# ---- extra tokens for 3.7
segs = []; prev = None
for e in Es:
    v = Bv(e, 40); m = int(v.argmax()) + 1
    if m != prev: segs.append([m, e, e]); prev = m
    else: segs[-1][2] = e
def fe(e): return '%g' % e
V['hiseq'] = ', '.join(('and ' if i == len(segs) - 1 and len(segs) > 1 else '') + ('Mixture~%d up to %s~MeV' % (m, fe(e1)) if i == 0 else 'Mixture~%d from %s~MeV' % (m, fe(e0)) if i == len(segs) - 1 else 'Mixture~%d between %s and %s~MeV' % (m, fe(e0), fe(e1))) for i, (m, e0, e1) in enumerate(segs))
V['ebf5_hi'] = '%.3g' % Bv(5.0, 40).max()
old = s
def para(start, end=r'\n\n'):
    a = old.index(start); b = old.index('\n\n', a); return old[a:b]
p_intro_ebf = para(r'The narrow-beam attenuation coefficients quantify removal from the primary beam')
p_limit = para(r'Two limitations bound the interpretation of these results.')
p_equalmfp = para(r'These comparisons are made at equal mfp rather than equal physical thickness.')
p_lowE_old = para(r'The low-energy side of the buildup curve is also physically relevant.')
p_compare_note = para(r'The present low-$Z$ mixtures should not be read as competitors')
p_regime = para(r'The rise of all three thickness scales toward high energy follows the interaction decomposition')
p_regime = p_regime.replace('these CHON-only systems necessarily', 'these low-$Z$ systems necessarily')
p_mfp_intro = para(r'Mean free path (MFP), half-value layer (HVL)'); p_rpe_intro = para(r'Radiation protection efficiency (RPE) expresses the same LAC information')
p_maclac_intro = para(r'Figure~\ref{fig:maclac} shows the full energy dependence of MAC and LAC')

sec33 = (r'''\subsection{Mass-Normalized Versus Density-Weighted Attenuation}\label{sec:maclac}

''' + p_maclac_intro + '\n\n' + figML + r'''

At 0.01~MeV the two effects act together. MAC spans @MAC001@~cm\textsuperscript{2}/g (@sp001@\%), while LAC spans @LAC001@~cm\textsuperscript{$-1$} (@lacsp001@\%). Mixture~3 has both the largest MAC and the largest LAC, and Mixture~2 the smallest of both. Mixtures~1--3 contain the same BSA fraction and therefore the same nitrogen and sulfur, so their differences at this energy arise from the balance between oxygen on one side and carbon and hydrogen on the other: Mixture~3 is the most oxygen-rich of the three and Mixture~2 the most carbon-rich, and the contrasts of Table~\ref{tbl:sensitivity} (oxygen raising and hydrogen lowering the MAC relative to carbon) account for a difference of the observed size. Mixture~2 also has the lowest density, which reinforces its low LAC.

The distinction strengthens beyond the photoelectric region. At 0.1~MeV, MAC values have converged to @MAC010@~cm\textsuperscript{2}/g, a spread of only @sp010@\%, while LAC still spans @LAC010@~cm\textsuperscript{$-1$}, a @lacsp010@\% difference. At 10~MeV the Geant4 MAC spread of @sp10@\% is comparable to the counting noise (the noise-free XCOM value is @xsp10@\%), whereas LAC spans @LAC1000@~cm\textsuperscript{$-1$}, a spread of @lacsp1000@\%. Once mass-normalized attenuation is nearly composition-independent, the residual separation in attenuation per unit length is therefore governed by density rather than by any large difference in intrinsic interaction probability per unit mass.

An equal-transmission comparison makes the practical consequence clear. At 0.1~MeV, Mixture~10, the densest, requires @thin@\% less physical thickness than Mixture~2, the least dense, to produce the same primary-beam transmission, because its LAC is larger. The required areal mass, however, scales as $1/\mathrm{MAC}$ at fixed transmission, and Mixture~10 has the slightly smaller MAC, so its areal mass is @am@\% $\pm$ @amu@\% higher (@amx@\% from the noise-free XCOM coefficients). Increased density therefore delivers a compactness advantage rather than a mass-efficiency advantage. The sign of the MAC difference is that expected from the hydrogen contrast of Section~\ref{sec:composition}: between 0.1 and 3~MeV each weight percent of hydrogen replacing carbon raises the MAC by @dHmin@--@dHmax@\%, and Mixture~2 carries the highest hydrogen mass fraction of the library (@H2@ against @Hm10@ for Mixture~10), which raises its electron density per unit mass ($Z/A\approx1$ for hydrogen against $\approx0.5$ for C, N, and O). The difference between these two mixtures is of the size that hydrogen alone predicts, but the library cannot by itself prove that hydrogen is the cause, since the elemental fractions vary together; the decorrelated analysis of Section~\ref{sec:composition} is what identifies hydrogen as the controlling element in this energy range.

''' + p_compare_note + r'''

''')
sec34 = (r'''\subsection{Thickness-Derived Attenuation Scales: MFP, HVL, and TVL}\label{sec:mfphvltvl}

''' + p_mfp_intro + '\n\n' + figMF + r'''

At 0.01~MeV, where photoelectric absorption dominates, HVL spans only @HVL001@~cm. By 0.1~MeV it has increased to @HVL010@~cm, reflecting the rapid loss of photoelectric absorption and the transition into Compton-dominated transport, and at 15~MeV it reaches @HVL1500@~cm (MFP @MFP1500@~cm, TVL @TVL1500@~cm). MFP and TVL reproduce the same evolution because their ratios to HVL are fixed by definition.

A second feature is more informative for the composition analysis. Once the MAC values have converged in the Compton region, the residual separation among the thickness curves persists because LAC retains the density factor. Mixture~2, the least dense, has the largest thickness scales at all @nE@ energies, and Mixture~10, the densest, the smallest at @n10@ of them (Mixture~3 at the other eight) --- the same compactness effect quantified in Section~\ref{sec:maclac}.

''' + p_regime + r'''

''')
sec35 = (r'''\subsection{Radiation Protection Efficiency at Equal Physical Thickness}\label{sec:rpe}

''' + p_rpe_intro + '\n\n' + figR + r'''

At 0.01~MeV the LAC range of @LAC001@~cm$^{-1}$ corresponds to an RPE of @RPE001@\% for a 1~cm slab. Efficiency then drops sharply as photoelectric absorption gives way to Compton scattering: at 0.1~MeV the LAC interval @LAC010@~cm$^{-1}$ gives an RPE of @RPE010@\%. By 1~MeV the corresponding RPE is @RPE100@\%, falling further to @RPE1000@\% at 10~MeV and @RPE1500@\% at 15~MeV. These values are the nonlinear equal-thickness representation of the attenuation coefficient discussed in Sections~\ref{sec:maclac} and~\ref{sec:mfphvltvl}.

The inter-mixture ordering remains physically informative. Mixture~10 has the highest RPE at @n10rpe@ of the @nE@ energies (Mixture~3 at the low-energy points and at 6 and 12~MeV) and Mixture~2 the lowest at all of them, because RPE is monotonic in LAC at fixed $x$. The magnitude of that separation is specific to the chosen thickness, since Eq.~\eqref{eq:rpe} is nonlinear in $\mu x$ and saturates at large optical depth.

''')
sec36 = (r'''\subsection{Energy-Dependent Effective Atomic Number}\label{sec:zeff}

The effective atomic number $Z_{\mathrm{eff}}(E)$ provides a compact interaction-weighted descriptor of the changing elemental response, and supports the process-resolved analysis of Section~\ref{sec:composition}. Figure~\ref{fig:zeff} shows that $Z_{\mathrm{eff}}$ is strongly energy dependent even though the elemental compositions are fixed.

At 0.01~MeV, $Z_{\mathrm{eff}}$ ranges from @Zlo001@ for Mixture~10 to @Zhi001@ for Mixture~4. The ordering follows sulfur (equivalently nitrogen and BSA) and not oxygen: the correlation of $Z_{\mathrm{eff}}$ with the sulfur fraction is $@rSZ001@$ at 0.01~MeV, against $@rOZ001@$ for oxygen (Table~\ref{tbl:correlation}), as expected because photoelectric absorption weights the highest-$Z$ constituent, sulfur, most heavily and Mixture~10 contains none. As the photon energy enters the Compton-dominated region, $Z_{\mathrm{eff}}$ decreases rapidly for every mixture and reaches @Z1lo@--@Z1hi@ near 1~MeV. The ranking changes: from 0.02~MeV upward Mixture~2 is the lowest at every energy, the oxygen coefficient becomes positive ($@rOZ003@$ at 0.03~MeV and $@rOZ100@$ at 1~MeV) and the sulfur coefficient negative ($@rSZ100@$ at 1~MeV), showing that no single composition variable fixes the $Z_{\mathrm{eff}}$ ordering across all interaction regimes. At high energy the growing pair-production contribution raises $Z_{\mathrm{eff}}$ again, to @Zlo1500@--@Zhi1500@ at 15~MeV.

''' + figZ + r'''

Table~\ref{tbl:zeff} lists the calculated values at selected energies. The absolute inter-mixture range remains modest: the spread is @Zsp001@ in $Z_{\mathrm{eff}}$ at 0.01~MeV, @Zsp100@ near 1~MeV, and @Zsp1500@ at 15~MeV. Mixture~4 has the highest $Z_{\mathrm{eff}}$ at every energy, and above 0.02~MeV Mixture~2 has the lowest; for example, at 10~MeV Mixture~4 exceeds Mixture~10 (@z10_4@ against @z10_10@). Composition tuning therefore changes the interaction-weighted atomic response measurably, particularly in the photoelectric region, but does not move the system into a qualitatively different atomic-number regime. The minimum near 1~MeV marks the energy at which the response is most purely Compton-like, and the subsequent rise identifies the onset of pair production as a second $Z$-sensitive mechanism.

Applying the same definition and cross-section data to water and PMMA places the library in context. Figure~\ref{fig:zeff_baseline} shows that at 0.01~MeV water ($Z_{\mathrm{eff}}=@Zwat001@$) lies above the complete range of the mixtures (@Zlo001@--@Zhi001@), whereas PMMA (@Zpm001@) lies below it, and PMMA remains below the range at every energy. Water lies within the range at 0.02~MeV (@Zwat020@ against @Zlo020@--@Zhi020@), coincides with its lower edge at 0.03~MeV (@Zwat030@), and lies below it from 0.04~MeV upward: at 1~MeV the mixtures span @Z1lo@--@Z1hi@ against @Zwat100@ for water and @Zpm100@ for PMMA, and at 15~MeV the corresponding values are @Zlo1500@--@Zhi1500@, @Zwat1500@, and @Zpm1500@. The change in ordering follows the interaction-weighted elemental balance. Water contains @Hwater@~wt\% hydrogen against @Hmean@~wt\% (@Hlo@--@Hhi@~wt\%) in the mixtures, and hydrogen ($Z=1$) lowers $Z_{\mathrm{eff}}$ most strongly once Compton scattering dominates, whereas the high oxygen content of water (88.8~wt\%) gives it the largest $Z_{\mathrm{eff}}$ at 0.01~MeV. The comparison therefore reinforces the central result that effective atomic number is energy dependent: from 0.04~MeV upward the mixtures as a group lie above both references, whereas at the lowest energies water lies above them and PMMA below them.

''' + figZB + '\n\n' + tabZ + r'''

''')
sec37 = (r'''\subsection{Exposure Buildup and the Primary--Secondary Photon Balance}\label{sec:ebf}

''' + p_intro_ebf + r'''

Figure~\ref{fig:ebf} shows the characteristic buildup maximum near 0.08--0.1~MeV, in the energy region where Compton scattering redistributes rather than absorbs photons. The peak grows strongly with optical depth, from @pk5@ at 5~mfp to @pk10@ at 10~mfp and to $@pk40lo@\times10^{4}$--$@pk40hi@\times10^{4}$ at 40~mfp. Large numerical EBF values at large $d$ do not imply large absolute transmission, because the denominator of the buildup ratio is itself exponentially small. At 0.1~MeV and 40~mfp, for instance, $e^{-40}\approx4.25\times10^{-18}$, so that multiplying by the calculated buildup factors gives total-response proxies of only about $@proxlo@\times10^{-13}$ for Mixture~@m010lo@ and $@proxhi@\times10^{-13}$ for Mixture~@m010hi@.

''' + figE + r'''

The low-energy side of the buildup curve is also physically relevant. At 40~mfp, EBF ranges from @E400015@ at 0.015~MeV and @E40002@ at 0.02~MeV to @E40003@ at 0.03~MeV before reaching the much larger maximum near 0.08--0.1~MeV. Photoelectric absorption dominates primary removal in this region, but photons that do scatter accumulate over many interaction lengths, so buildup begins to rise before the Compton-dominated plateau in the total MAC is reached.

Mixture separation is largest neither where the absolute EBF is largest nor in a fixed order. At 0.04~MeV and 40~mfp the relative spread reaches @ebfsp004@\%, from @b004lo@ for Mixture~@m004lo@ to @b004hi@ for Mixture~@m004hi@. At 0.1~MeV and 40~mfp, Mixture~@m010hi@ has an EBF of @b010hi@ and Mixture~@m010lo@ an EBF of @b010lo@, a @ebfsp010@\% spread corresponding to a @red010@\% reduction for Mixture~@m010lo@. The mixture with the highest EBF at 40~mfp is @hiseq@. Between 0.04 and about 1~MeV the highest buildup therefore belongs to Mixture~10, the formulation without sulfur or nitrogen, and the lowest to Mixture~4. The G--P ranking follows the equivalent atomic number rather than the direct $Z_{\mathrm{eff}}$: at 0.1~MeV the EBF at 40~mfp decreases monotonically with $Z_{\mathrm{eq}}$ across all ten mixtures, with $Z_{\mathrm{eq}}=@Zq10@$ for Mixture~10 and @Zq4@ for Mixture~4, and sulfur and nitrogen raise $Z_{\mathrm{eq}}$; by contrast, Mixture~10 has the third-highest direct $Z_{\mathrm{eff}}$ at that energy. The direction of this result matters for shielding: the mixture with the largest LAC also produces the most secondary buildup, so its advantage in primary attenuation is partly offset in the total response. Above 1~MeV the differences become small (@sp1e@\% at 1~MeV, @sp10e@\% at 10~MeV, and @sp15e@\% at 15~MeV, where Mixture~@m15hi@ has an EBF of @b15hi@ against @b15lo@ for Mixture~@m15lo@), and the buildup ranking changes with energy, as do the $Z_{\mathrm{eff}}$ and high-energy MAC rankings.

''' + p_equalmfp.replace('The present data therefore establish that secondary-photon buildup does not follow a fixed composition ranking, while determining which mixture minimises total broad-beam response at a specified physical thickness would require evaluating', 'The present data therefore establish that secondary-photon buildup does not follow a fixed composition ranking, while determining which mixture minimises total broad-beam response at a specified physical thickness requires evaluating') + r'''

''' + figBB + r'''

Figure~\ref{fig:broadbeam} makes this comparison for two fixed conditions. At equal thickness (10~cm) the best-to-worst spread in primary transmission is @bbpT01@\% at 0.1~MeV, but it falls to @bbtT01@\% in the total response, because the mixture that attenuates the primary beam most strongly also builds up most; over 0.1--1~MeV the primary spread exceeds the total spread by a factor of @bbrTlo@--@bbrThi@, and above 1.3~MeV the total spread is at most @bbhiT@\%. The mixture that transmits least of the primary beam differs from the one that transmits least of the total response at @bbdiffT@ of the @bbnT@ energies for which the G--P range is valid (0.015--@bbEmaxT@~MeV): at 0.1~MeV it is Mixture~@bbT_bestp@ for the primary beam and Mixture~@bbT_bestb@ for the total response. At equal areal mass (20~g/cm\textsuperscript{2}) the ranking by primary attenuation follows the MAC and the spread is only @bbpM01@\% at 0.1~MeV, smaller than that of the total response (@bbtM01@\%), which is set mainly by buildup. The two rankings are essentially unrelated (median Spearman coefficient @bbspMmed@ over the @bbnM@ energies, minimum $@bbspMmin@$), and the best mixture differs at @bbdiffM@ of the @bbnM@ energies. Above 1.3~MeV the primary spread is at most @bbhipM@\% and that of the total response at most @bbhiM@\%.

''' + p_limit + r'''

''')
new_results = sec32 + sec33 + sec34 + sec35 + sec36 + sec37
a = s.index(r'\subsection{Composition Sensitivity Across Photon-Interaction Regimes}')
if a < 0: raise SystemExit('start not found')
end_marker = s.index(r'\section{Conclusions}')
b = s.rfind('\n% ====', 0, end_marker)
s = s[:a] + new_results + s[b:]
for k, v in V.items(): s = s.replace('@%s@' % k, str(v))
left = re.findall(r'@\w+@', s); assert not left, left[:10]
open(TEX, 'w').write(s); print('results rewritten; chars', len(s))
