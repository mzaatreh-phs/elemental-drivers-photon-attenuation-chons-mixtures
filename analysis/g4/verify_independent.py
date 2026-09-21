#!/usr/bin/env python3
"""Independent re-derivation (pure Python, no numpy, no shared code with g4analysis.py) of pooled MAC, uncertainty, XCOM deviation
and the MAC spread at each complete energy, straight from the raw per-run CSVs. Compares with analysis/g4/results/mac_pooled.csv.
Exit code 1 if any value differs by more than the tolerance."""
import csv, glob, os, math, sys
R = os.path.expanduser('~/CHON/recalc_2026-09-19')
raw = {}
for d in glob.glob(f'{R}/runs/points/s[12]_Mix_*_E*'):
    if not os.path.exists(d + '/done.flag'): continue
    r = list(csv.DictReader(open(glob.glob(d + '/MAC_results_*.csv')[0])))[0]
    key = (int(r['material_id'].split('_')[1]), round(float(r['Energy(MeV)']), 6))
    mac_txt = open(f"{R}/macros/points/{os.path.basename(d)}.mac").read()          # density exactly as passed to Geant4
    rho_macro = float(mac_txt.split('/det/setDensity')[1].split()[0])
    raw.setdefault(key, []).append((int(r['n_incident']), int(r['n_transmitted_strict']), rho_macro, float(r['Thikness(cm)'])))
mine = {}
for (m, E), v in raw.items():
    if len(v) != 2: continue
    N = sum(a[0] for a in v); k = sum(a[1] for a in v); rho = v[0][2]; x = v[0][3]
    T = k / N; mac = -math.log(T) / (rho * x); u = math.sqrt((1 - T) / (N * T)) / (rho * x)
    mine[(m, E)] = (mac, u, 100 * (1 - T))
worst = {'mac': 0, 'u': 0, 'rpe': 0}; n = 0
for r in csv.DictReader(open(f'{R}/analysis/g4/results/mac_pooled.csv')):
    key = (int(r['mixture']), round(float(r['E_MeV']), 6))
    if key not in mine: print('missing in independent set:', key); sys.exit(1)
    mac, u, rpe = mine[key]
    # allowed difference = 6e-7 (CSV prints 6 decimals) + 8e-7*MAC (macro density printed to 6 decimals, relative <= 5e-7/1.26)
    worst['mac'] = max(worst['mac'], abs(mac - float(r['MAC_G4'])) / (6e-7 + 8e-7 * mac))
    worst['u'] = max(worst['u'], abs(u - float(r['u_MAC'])))
    worst['rpe'] = max(worst['rpe'], abs(rpe - float(r['RPE_1cm_pct'])))
    n += 1
# spread check per energy
sp_bad = 0
import json
t4 = {round(x['E'], 6): x for x in json.load(open(f'{R}/analysis/g4/results/table4_spreads.json'))}
for E, x in t4.items():
    v = [mine[(m, E)][0] for m in range(1, 11)]
    sp = 100 * (max(v) / min(v) - 1)
    if abs(sp - x['MAC']['spread_pct']) > 1e-3: sp_bad += 1
print(f'compared {n} (mixture, energy) points; worst MAC diff / allowed {worst["mac"]:.2f}, u(MAC) {worst["u"]:.1e}, RPE {worst["rpe"]:.1e} pct-pts; spread mismatches: {sp_bad}')
# tolerances = half a unit in the last printed digit (6 decimals for MAC/u, 4 for RPE) plus a little margin
sys.exit(0 if worst['mac'] < 1.0 and worst['u'] < 6e-7 and worst['rpe'] < 6e-5 and sp_bad == 0 else 1)
