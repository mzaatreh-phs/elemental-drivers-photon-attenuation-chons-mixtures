"""Core library for the CHON Geant4 recalculation analysis (fixed 1.0 cm slab, 1e6 photons/run, seeds s1/s2).

Everything here is derived from raw counts (n_incident, n_transmitted_strict) so uncertainties are binomial and traceable.
Run with:  ~/CHON/recalc_2026-09-19/.venv/bin/python analysis/g4/analyze.py
"""
import csv, json, glob, os, math
import numpy as np

R = os.path.expanduser('~/CHON/recalc_2026-09-19')
RUNS = f'{R}/runs/points'
XC = f'{R}/analysis/xcom'
OUT = f'{R}/analysis/g4'
MIXES = list(range(1, 11))

# ---------------------------------------------------------------- inputs
def load_inputs():
    d = json.load(open(f'{R}/inputs/inputs_exact.json'))
    mix = d['mixtures']
    comp = {int(m): {e: x['composition'].get(e, 0.0) for e in 'CHONS'} for m, x in mix.items()}
    rho = {int(m): x['density'] for m, x in mix.items()}
    bsa = {int(m): x['weights']['BSA'] for m, x in mix.items()}
    grid = sorted(d['production_energies'] + d['validation_energies'])
    return comp, rho, bsa, grid, sorted(d['validation_energies'])

COMP, RHO, BSA, GRID, VALID = load_inputs()

# ---------------------------------------------------------------- raw runs
def load_runs():
    """Return {(seed_label, mix, E): dict(n_inc, n_tr, x, rho)} for finished jobs; also list of problems."""
    runs, problems = {}, []
    for d in sorted(glob.glob(f'{RUNS}/s[12]_Mix_*_E*')):
        if not os.path.exists(f'{d}/done.flag'):
            continue
        f = glob.glob(f'{d}/MAC_results_*.csv')
        rows = list(csv.DictReader(open(f[0]))) if f else []
        if len(rows) != 1:
            problems.append((os.path.basename(d), 'rows=%d' % len(rows))); continue
        r = rows[0]
        lab = os.path.basename(d).split('_')[0]
        m = int(r['material_id'].split('_')[1]); E = float(r['Energy(MeV)'])
        rho = float(r['density_g_cm3'])
        if abs(rho - RHO[m]) > 2e-5:
            problems.append((os.path.basename(d), f'density {rho} != {RHO[m]}'))
        runs[(lab, m, round(E, 6))] = dict(n_inc=int(r['n_incident']), n_tr=int(r['n_transmitted_strict']),
                                          x=float(r['Thikness(cm)']), rho=RHO[m], csv_mac=float(r['MAC_strict(cm2/g)']))
    return runs, problems

def mac_from_counts(n_inc, n_tr, x, rho):
    """MAC (cm2/g), its binomial standard uncertainty and T."""
    n_inc = np.asarray(n_inc, float); n_tr = np.asarray(n_tr, float)
    T = n_tr / n_inc
    mac = -np.log(T) / (rho * x)
    u = np.sqrt((1 - T) / (n_inc * T)) / (rho * x)
    return mac, u, T

def complete_energies(runs):
    """Energies for which all 10 mixtures x 2 seeds are finished."""
    have = {}
    for (lab, m, E) in runs:
        have.setdefault(E, set()).add((lab, m))
    need = {(l, m) for l in ('s1', 's2') for m in MIXES}
    return sorted(E for E, s in have.items() if s >= need)

def build_tables(runs, energies):
    """Arrays indexed [mix, energy]: per-seed and pooled MAC/u/T."""
    nm, ne = len(MIXES), len(energies)
    A = {k: np.zeros((nm, ne)) for k in ['n1', 'k1', 'n2', 'k2']}
    x = np.zeros((nm, ne)); rho = np.zeros((nm, ne))
    for i, m in enumerate(MIXES):
        for j, E in enumerate(energies):
            a = runs[('s1', m, E)]; b = runs[('s2', m, E)]
            A['n1'][i, j], A['k1'][i, j], A['n2'][i, j], A['k2'][i, j] = a['n_inc'], a['n_tr'], b['n_inc'], b['n_tr']
            assert a['x'] == b['x']; x[i, j] = a['x']; rho[i, j] = RHO[m]
    D = {}
    D['E'] = np.array(energies); D['x'] = x; D['rho'] = rho
    D['mac1'], D['u1'], D['T1'] = mac_from_counts(A['n1'], A['k1'], x, rho)
    D['mac2'], D['u2'], D['T2'] = mac_from_counts(A['n2'], A['k2'], x, rho)
    D['nP'] = A['n1'] + A['n2']; D['kP'] = A['k1'] + A['k2']
    D['mac'], D['u'], D['T'] = mac_from_counts(D['nP'], D['kP'], x, rho)
    D['lac'] = D['mac'] * rho; D['u_lac'] = D['u'] * rho
    D['sigT'] = np.sqrt(D['T'] * (1 - D['T']) / D['nP'])          # binomial sd of pooled T
    return D

# ---------------------------------------------------------------- reference data (XCOM / EpiXS-equivalent)
def load_xcom_mac():
    out = {}
    for r in csv.DictReader(open(f'{XC}/process_fractions.csv')):
        if r['set'] == 'new':
            out[(int(r['material'][3:]), round(float(r['E_MeV']), 6))] = float(r['mac_total_cm2_g'])
    return out

def load_zeff():
    z, ref = {}, {}
    for r in csv.DictReader(open(f'{XC}/zeff.csv')):
        E = round(float(r['E_MeV']), 6)
        if r['set'] == 'new': z[(int(r['material'][3:]), E)] = float(r['Zeff_direct'])
        elif r['set'] == 'ref': ref[(r['material'], E)] = float(r['Zeff_direct'])
    return z, ref

# ---------------------------------------------------------------- statistics
def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a - a.mean(); b = b - b.mean()
    den = math.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else float('nan')

def spread_pct(v):
    v = np.asarray(v, float)
    return 100 * (v.max() / v.min() - 1)

def mc_draw(mean, sd, n=4000, seed=12345):
    """Independent normal draws around per-mixture means -> array [n, mixtures]."""
    rng = np.random.default_rng(seed)
    return mean[None, :] + sd[None, :] * rng.standard_normal((n, len(mean)))

def interval(x, lo=16, hi=84):
    return tuple(np.percentile(x, [lo, 50, hi]))
