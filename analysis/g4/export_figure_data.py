#!/usr/bin/env python3
"""Export the numbers behind every figure, ONE FOLDER PER FIGURE (folder name = PNG name in figs_new/).

Inside a folder every CSV is one panel (or one part of a panel: curves / error bars), named
    panelA_<what>.csv, panelB_<what>.csv ...   (panel order = left-to-right, top-to-bottom in the PNG)
Wide files: first column = x (Energy_MeV), then one column per curve (Mix1..Mix10 or the named series).
Each folder has a README.txt (which PNG, which panel, x column, y columns, axis scales); figure_data/README.txt is the index.
Also writes figure_data_all_sheets.xlsx (INDEX sheet + one sheet per CSV).
"""
import os, sys, csv, json, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from g4analysis import *
DEST = os.path.expanduser('~/CHON/CHOM_submitted version/figure_data/')
SENS = R + '/analysis/sensitivity/'; BB = R + '/analysis/broadbeam/'
# start clean: remove previously exported figure folders / old flat files (all regenerated below)
for f in os.listdir(DEST) if os.path.isdir(DEST) else []:
    p = DEST + f
    if os.path.isdir(p) and f in ('fig_elemental_contrast', 'fig_validation_residuals'): continue      # written by origin_style/make_new_figures.py
    if os.path.isdir(p): shutil.rmtree(p)
    elif f.endswith('.csv') or f.endswith('.xlsx') or f == 'README.txt': os.remove(p)
os.makedirs(DEST, exist_ok=True)

runs, _ = load_runs(); energies = complete_energies(runs); assert len(energies) == 38, len(energies)
D = build_tables(runs, energies); E = D['E']; xc = load_xcom_mac(); zeffd, zref = load_zeff()
Xc = np.array([[xc[(m, round(e, 6))] for e in E] for m in MIXES])
HDR = ['Energy_MeV'] + [f'Mix{m}' for m in MIXES]
INDEX = []      # (figure folder, file, description)
SHEETS = {}     # sheet name -> rows
ALIAS = {'fig_mac_lac': 'mac_lac', 'fig_mfp_hvl_tvl': 'mfp_hvl_tvl', 'fig_rpe': 'rpe', 'fig_validation_deviation': 'validation', 'fig_seed_check': 'seed_check',
         'fig_zeff': 'zeff', 'fig_zeff_baseline': 'zeff_base', 'fig_process_fractions': 'process', 'fig_ebf': 'ebf', 'fig_composition_sensitivity': 'comp_sens',
         'fig_elemental_sensitivity': 'elem_sens', 'fig_library_vs_ensemble_correlation': 'lib_vs_ens', 'fig_element_importance': 'elem_import',
         'fig_exchange_confirmation': 'exchange', 'fig_broadbeam_primary_vs_total': 'bb_prim_vs_tot', 'fig_broadbeam_spread': 'bb_spread', '.': 'Table2'}
def _sheet(folder, name):   # Excel sheet name = <figure alias>_<file name without 'panel'>, <= 31 chars
    nm = name.replace('panel', '').replace('_percent', '').replace('_relative_to_mean', '_rel')
    for pre in (ALIAS[folder] + '_', 'Table2_', 'validation_'): nm = nm[len(pre):] if nm.startswith(pre) else nm
    return f"{ALIAS[folder]}_{nm}"[:31]
def _out(folder, name, header, rows, what):
    os.makedirs(DEST + folder, exist_ok=True)
    with open(f'{DEST}{folder}/{name}.csv', 'w', newline='') as f: csv.writer(f).writerows([header] + rows)
    INDEX.append((folder, name + '.csv', what))
    sh = _sheet(folder, name); assert sh not in SHEETS, sh
    SHEETS[sh] = [header] + [[float(x) if _isnum(x) else x for x in r] for r in rows]
def _isnum(x):
    try: float(x); return not isinstance(x, str) or x.strip() != ''
    except (TypeError, ValueError): return False
def wide(folder, name, arr, what, energies_=None, fmt='%.8g'):
    en = E if energies_ is None else energies_
    _out(folder, name, HDR, [[f'{e:g}'] + [fmt % arr[i][j] for i in range(len(MIXES))] for j, e in enumerate(en)], what)
def table(folder, name, header, rows, what): _out(folder, name, header, rows, what)
def readme(folder, text):
    open(f'{DEST}{folder}/README.txt', 'w').write(text.strip() + '\n')
MIXNOTE = 'Columns Mix1..Mix10 = the ten mixtures (composition/density in Table2_mixture_composition_and_density.csv).'

# ---------------------------------------------------------------- Geant4 figures
F = 'fig_mac_lac'
wide(F, 'panelA_MAC_curves', D['mac'], 'MAC (Geant4) - the curves of the left panel')
wide(F, 'panelA_MAC_errorbars_1sigma', D['u'], 'MAC 1-sigma binomial uncertainty (error bars of left panel; smaller than the markers in the PNG)')
wide(F, 'panelA_MAC_XCOM_reference', Xc, 'XCOM total MAC with coherent scattering, same energies - NOT plotted in the PNG, optional overlay')
wide(F, 'panelB_LAC_curves', D['lac'], 'LAC = MAC x density (Geant4) - the curves of the right panel')
wide(F, 'panelB_LAC_errorbars_1sigma', D['u_lac'], 'LAC 1-sigma uncertainty (error bars of right panel)')
readme(F, f"""figs_new/fig_mac_lac.png - Geant4, 38 energies 0.01-15 MeV, both seeds pooled. Log-log axes.
PANEL A (left):  x = Energy_MeV,  y = MAC (cm2/g)  -> panelA_MAC_curves.csv (10 curves). Error bars: panelA_MAC_errorbars_1sigma.csv.
PANEL B (right): x = Energy_MeV,  y = LAC (1/cm)   -> panelB_LAC_curves.csv (10 curves). Error bars: panelB_LAC_errorbars_1sigma.csv.
{MIXNOTE}""")

F = 'fig_mfp_hvl_tvl'
wide(F, 'panelA_MFP_cm', 1 / D['lac'], 'mean free path = 1/LAC (Geant4)')
wide(F, 'panelB_HVL_cm', np.log(2) / D['lac'], 'half-value layer = ln2/LAC (Geant4)')
wide(F, 'panelC_TVL_cm', np.log(10) / D['lac'], 'tenth-value layer = ln10/LAC (Geant4)')
readme(F, f"""figs_new/fig_mfp_hvl_tvl.png - Geant4. Three panels, log-log axes, x = Energy_MeV.
PANEL A (left):   y = MFP (cm)  -> panelA_MFP_cm.csv
PANEL B (middle): y = HVL (cm)  -> panelB_HVL_cm.csv
PANEL C (right):  y = TVL (cm)  -> panelC_TVL_cm.csv
{MIXNOTE}""")

F = 'fig_rpe'
wide(F, 'curves_attenuation_1cm_percent', 100 * (1 - D['T']), 'primary-beam attenuation (1 - T) x 100 for the 1.0 cm slab (Geant4) - the curves')
wide(F, 'errorbars_1sigma_percent', 100 * D['sigT'], '1-sigma of the attenuation (error bars)')
readme(F, f"""figs_new/fig_rpe.png - Geant4, 1.0 cm slab, single panel, x log axis (y linear).
x = Energy_MeV, y = primary-beam attenuation (%) -> curves_attenuation_1cm_percent.csv. Error bars: errorbars_1sigma_percent.csv.
{MIXNOTE}""")

F = 'fig_validation_deviation'
wide(F, 'points_deviation_percent', 100 * (D['mac'] / Xc - 1), 'signed deviation 100 x (MAC_Geant4/MAC_XCOM - 1) - the plotted points')
wide(F, 'errorbars_1sigma_percent', 100 * D['u'] / Xc, 'error bar length (1 sigma) of each point')
table(F, 'validation_energies_MeV', ['Energy_MeV'], [[f'{v:g}'] for v in VALID], 'the six validation energies (dotted vertical lines in the PNG)')
readme(F, f"""figs_new/fig_validation_deviation.png - Geant4 vs XCOM, single panel, x log axis.
x = Energy_MeV, y = deviation (%) -> points_deviation_percent.csv (10 series), error bars: errorbars_1sigma_percent.csv.
Horizontal line at y = 0. Dotted vertical lines at validation_energies_MeV.csv.
In the PNG the ten series are shifted horizontally by at most +-1.2% for visibility (markers would overlap); the data are at the exact energies.
{MIXNOTE}""")

z = (D['mac1'] - D['mac2']) / np.sqrt(D['u1'] ** 2 + D['u2'] ** 2)
F = 'fig_seed_check'
wide(F, 'z_values_all_points', z, 'z = (MAC_seed1 - MAC_seed2)/sqrt(u1^2+u2^2) for every mixture and energy - raw data of the histogram')
cnt, edges = np.histogram(z.ravel(), bins=np.linspace(-4, 4, 33), density=True)
table(F, 'histogram_and_N01_curve', ['bin_left', 'bin_right', 'bin_center', 'density', 'N01_density'],
      [[edges[i], edges[i + 1], (edges[i] + edges[i + 1]) / 2, cnt[i], np.exp(-((edges[i] + edges[i + 1]) / 2) ** 2 / 2) / np.sqrt(2 * np.pi)] for i in range(len(cnt))],
      'histogram bins (bar height = density) and the standard-normal curve (black line in the PNG)')
readme(F, """figs_new/fig_seed_check.png - single panel histogram, x linear.
Bars: histogram_and_N01_curve.csv, x = bin_center (width bin_right-bin_left), y = density (32 bins from -4 to 4, all 380 z values).
Line: same file, x = bin_center, y = N01_density (standard normal N(0,1)).
Raw z values behind the histogram: z_values_all_points.csv (rows = energy, columns = mixture).""")

# ---------------------------------------------------------------- XCOM / G-P figures
F = 'fig_zeff'
Zm = np.array([[zeffd[(m, round(e, 6))] for e in E] for m in MIXES]); wide(F, 'curves_Zeff', Zm, 'Zeff of each mixture (XCOM direct method)')
readme(F, f"""figs_new/fig_zeff.png - XCOM direct method, exact CHONS compositions, single panel, x log axis.
x = Energy_MeV, y = Zeff -> curves_Zeff.csv (10 curves).
{MIXNOTE}""")
F = 'fig_zeff_baseline'
lo, hi = Zm.min(0), Zm.max(0)
table(F, 'band_and_reference_Zeff', ['Energy_MeV', 'CHONS_band_min', 'CHONS_band_max', 'Water', 'PMMA'],
      [[f'{e:g}', '%.6g' % lo[j], '%.6g' % hi[j], '%.6g' % zref[('water', round(e, 6))], '%.6g' % zref[('PMMA', round(e, 6))]] for j, e in enumerate(E)],
      'shaded band (min..max Zeff over the ten mixtures) + water and PMMA lines')
readme(F, """figs_new/fig_zeff_baseline.png - single panel, x log axis.
x = Energy_MeV. Shaded band: fill between CHONS_band_min and CHONS_band_max. Lines: Water, PMMA (all in band_and_reference_Zeff.csv).""")

F = 'fig_process_fractions'
P = [r for r in csv.DictReader(open(f'{XC}/process_fractions.csv')) if r['set'] == 'new' and r['material'] == 'Mix1']
table(F, 'Mix1_process_fractions_percent', ['Energy_MeV', 'Photoelectric', 'Incoherent', 'Coherent', 'Pair_production'],
      [[f"{float(r['E_MeV']):g}"] + ['%.6g' % (100 * float(r[k])) for k in ('F_photoelectric', 'F_incoherent', 'F_coherent', 'F_pair')] for r in P],
      'Mixture 1: percentage of the total MAC from each process (4 curves)')
table(F, 'crossover_keV_all_mixtures', ['Mixture', 'Photoelectric_incoherent_crossover_keV'],
      [[r['material'], r['crossover_keV']] for r in csv.DictReader(open(f'{XC}/crossover.csv')) if r['set'] == 'new'],
      'photoelectric = incoherent energy per mixture; the dashed vertical line in the PNG is the Mix1 value')
readme(F, """figs_new/fig_process_fractions.png - XCOM, Mixture 1 only, single panel, x log axis.
x = Energy_MeV, y = fraction of total MAC (%) -> Mix1_process_fractions_percent.csv (4 curves: photoelectric, incoherent, coherent, pair production).
Dashed vertical line at the Mix1 crossover, in keV (divide by 1000 for MeV) -> crossover_keV_all_mixtures.csv (all ten mixtures listed).""")

F = 'fig_ebf'
EB = {}
for r in csv.DictReader(open(f'{XC}/exposure_buildup.csv')):
    if r['set'] == 'new': EB[(int(r['material'][3:]), round(float(r['E_MeV']), 6))] = r
E37 = [e for e in GRID if e >= 0.015]
for k, d in enumerate((5, 10, 15, 20, 25, 30, 35, 40)):
    arr = np.array([[float(EB[(m, round(e, 6))][f'B_{d}mfp']) for e in E37] for m in MIXES])
    wide(F, f'panel{k + 1}_{d}mfp_EBF', arr, f'exposure buildup factor at {d} mfp (G-P), panel {k + 1} of 8', E37, '%.6g')
readme(F, f"""figs_new/fig_ebf.png - G-P exposure buildup factor, 37 energies 0.015-15 MeV (G-P tables start at 0.015 MeV). 2 rows x 4 columns, log-log axes.
Panels in reading order: 1=5 mfp, 2=10, 3=15, 4=20 (top row); 5=25, 6=30, 7=35, 8=40 mfp (bottom row).
x = Energy_MeV, y = EBF -> panel<k>_<depth>mfp_EBF.csv (10 curves each).
{MIXNOTE}""")

# ---------------------------------------------------------------- composition sensitivity (Geant4 + XCOM)
F = 'fig_composition_sensitivity'
rowsA, rowsB = [], []
PR = {'Oxygen': [COMP[m]['O'] for m in MIXES], 'Sulfur': [COMP[m]['S'] for m in MIXES], 'Hydrogen': [COMP[m]['H'] for m in MIXES]}
for j, e in enumerate(E):
    s = np.array([spread_pct(d) for d in mc_draw(D['mac'][:, j], D['u'][:, j], n=1500)]); l, m_, h = interval(s)
    rowsA.append([f'{e:g}', '%.6g' % spread_pct(D['mac'][:, j]), '%.6g' % l, '%.6g' % h, '%.6g' % spread_pct(Xc[:, j])])
    rowsB.append([f'{e:g}'] + ['%.6g' % pearson(pv, D['mac'][:, j]) for pv in PR.values()] + ['%.6g' % pearson(pv, Xc[:, j]) for pv in PR.values()])
table(F, 'panelA_MAC_spread_percent', ['Energy_MeV', 'Geant4_spread', 'Geant4_noise_band_16pct', 'Geant4_noise_band_84pct', 'XCOM_spread'], rowsA,
      'inter-mixture MAC spread; line+markers = Geant4_spread, shaded band = 16..84 pct columns, dashed black = XCOM_spread')
table(F, 'panelB_Pearson_r_with_MAC', ['Energy_MeV', 'Oxygen_Geant4', 'Sulfur_Geant4', 'Hydrogen_Geant4', 'Oxygen_XCOM', 'Sulfur_XCOM', 'Hydrogen_XCOM'], rowsB,
      'Pearson r of MAC with O, S, H mass fraction; solid = Geant4 (markers), dashed = XCOM')
readme(F, """figs_new/fig_composition_sensitivity.png - two panels, x log axis (Energy_MeV).
PANEL A (left, log y): MAC inter-mixture spread (%) -> panelA_MAC_spread_percent.csv.
   Solid line with markers = Geant4_spread; shaded band between Geant4_noise_band_16pct and _84pct (68% Monte Carlo noise band); dashed black = XCOM_spread.
PANEL B (right, y from -1.05 to 1.05): Pearson r between MAC and the mass fraction of the element -> panelB_Pearson_r_with_MAC.csv.
   Solid + markers = *_Geant4 columns, dashed same colour = *_XCOM columns (Oxygen, Sulfur, Hydrogen).
Note: sulfur, nitrogen and BSA fraction are perfectly collinear in this ten-mixture library, so the sulfur curve is also the nitrogen and BSA curve.
Above 1 MeV the Geant4 spread is noise-limited; use the XCOM columns there.""")

# ---------------------------------------------------------------- novelty items 1-2 (sources copied from analysis/sensitivity and analysis/broadbeam)
def rd(fn):
    r = list(csv.reader(open(fn))); return r[0], r[1:]
h, a = rd(SENS + 'exchange_sensitivity_and_contribution.csv')
F = 'fig_elemental_sensitivity'
table(F, 'panelA_MAC_change_per_1wtpct_replacing_C_percent', ['Energy_MeV', 'H', 'N', 'O', 'S'], [[r[0]] + r[1:5] for r in a],
      'percent change of MAC per +1 wt% of the element replacing carbon (XCOM mixture rule, library-mean composition)')
table(F, 'panelB_share_of_MAC_percent_library_mean', ['Energy_MeV', 'C', 'H', 'N', 'O', 'S'], [[r[0]] + r[5:10] for r in a],
      'share of the total MAC contributed by each element (%), library-mean composition')
readme(F, f"""figs_new/fig_elemental_sensitivity.png - XCOM mixture rule, library-mean composition. Two panels, x log axis (Energy_MeV).
PANEL A (left, symlog y, linear within +-0.3): y = % change of MAC per +1 wt% of the element replacing C -> panelA_MAC_change_per_1wtpct_replacing_C_percent.csv (curves H, N, O, S). Horizontal line at 0.
PANEL B (right, log y): y = share of MAC (%) -> panelB_share_of_MAC_percent_library_mean.csv (curves C, H, N, O, S).
Source table: analysis/sensitivity/exchange_sensitivity_and_contribution.csv (headers there: {', '.join(h[:5])}, ...).""")

h, c = rd(SENS + 'ensemble_vs_library_correlation_and_importance.csv'); ix = {n: i for i, n in enumerate(h)}
F = 'fig_library_vs_ensemble_correlation'
for pn, el, nm in (('A', 'O', 'oxygen'), ('B', 'S', 'sulfur'), ('C', 'H', 'hydrogen')):
    table(F, f'panel{pn}_{nm}_r_with_MAC', ['Energy_MeV', 'r_ten_mixture_library', 'r_decorrelated_ensemble'], [[r[0], r[ix[f'r_library_{el}']], r[ix[f'r_ensemble_{el}']]] for r in c],
          f'Pearson r of MAC with {nm} mass fraction: ten-mixture library (red circles) vs decorrelated ensemble (blue squares)')
readme(F, """figs_new/fig_library_vs_ensemble_correlation.png - three panels, x log axis (Energy_MeV), y = Pearson r (fixed range -1.05..1.05), horizontal line at 0.
PANEL A (left) oxygen, PANEL B (middle) sulfur, PANEL C (right) hydrogen.
Each file: r_ten_mixture_library (red, circles) and r_decorrelated_ensemble (blue, squares).
Ensemble = 4000 random CHONS compositions inside the library's ranges (XCOM mixture rule).""")

F = 'fig_element_importance'
table(F, 'importance_percent', ['Energy_MeV', 'H', 'N', 'O', 'S'], [[r[0]] + [r[ix[f'importance_pct_{el}']] for el in 'HNOS'] for r in c],
      'importance of each element for MAC variability (%) in the decorrelated ensemble (standardized-beta based)')
readme(F, """figs_new/fig_element_importance.png - single panel, x log axis.
x = Energy_MeV, y = importance for MAC variability (%) -> importance_percent.csv (curves H, N, O, S).
CAUTION: only partly robust (see range test in extra_no_figure/range_test_importance.csv): hydrogen dominates 0.1-3 MeV in all scenarios; below ~0.06 MeV sulfur dominance
depends on the assumed element ranges; above ~8 MeV no element dominates robustly.""")

F = 'fig_exchange_confirmation'
ex = [r for r in csv.DictReader(open(SENS + 'exchange_confirmation.csv')) if r['variant'].startswith('d')]
VAR = [('dO', '+5 wt% O'), ('dN', '+5 wt% N'), ('dS', '+1.5 wt% S'), ('dH', '+1.5 wt% H'), ('dALL', 'all four')]
table(F, 'panelA_G4_vs_XCOM_MAC_change', ['variant', 'variant_label', 'Energy_MeV', 'XCOM_mixture_rule_change_percent_(x)', 'Geant4_change_percent_(y)', 'Geant4_1sigma_percent_(y_error)'],
      [[v, lab, r['E_MeV'], r['dMAC_XCOM_mixture_rule_pct'], r['dMAC_Geant4_pct'], r['u_pct']] for v, lab in VAR for r in ex if r['variant'] == v],
      'scatter of Geant4 vs XCOM change in MAC caused by each composition exchange (long format: 5 variants x 8 energies)')
Ev = sorted({float(r['E_MeV']) for r in ex})
zw = {(r['variant'], float(r['E_MeV'])): r['z'] for r in ex}
table(F, 'panelB_z_vs_energy', ['Energy_MeV'] + [lab for _, lab in VAR], [[f'{e:g}'] + [zw[(v, e)] for v, _ in VAR] for e in Ev],
      'z = (Geant4 - XCOM)/uncertainty at each energy, one column per variant (wide)')
readme(F, """figs_new/fig_exchange_confirmation.png - Geant4 confirmation of the XCOM mixture-rule prediction. Base = library-mean composition, 8 energies, 2 seeds.
(No plotting script was kept for this figure; these files are built from analysis/sensitivity/exchange_confirmation.csv.)
PANEL A (left, linear axes): x = XCOM_mixture_rule_change_percent_(x), y = Geant4_change_percent_(y), error bar = Geant4_1sigma_percent_(y_error), colour by variant
   -> panelA_G4_vs_XCOM_MAC_change.csv (long format; filter by column 'variant'). Dashed black line y = x.
PANEL B (right, x log axis): x = Energy_MeV, y = z, one series per variant -> panelB_z_vs_energy.csv (wide). Horizontal line at 0, dotted at +-2.
Variants: dO = +5 wt% O replacing C, dN = +5 wt% N, dS = +1.5 wt% S, dH = +1.5 wt% H, dALL = all four together.""")

# item 2: equal-thickness / equal-areal-mass broad-beam response
def wide_np(fn):
    hh, aa = rd(fn); return [r[0] for r in aa], np.array([[float(x) if x != '' else np.nan for x in r[1:]] for r in aa])
F = 'fig_broadbeam_primary_vs_total'
panel = 0
for row, (kind, lab) in enumerate((('primary', 'primary only'), ('broad', 'with buildup'))):
    for col, (tag, ttl, nm) in enumerate((('t10cm', 'equal thickness 10 cm', 'equal10cm'), ('m20gcm2', 'equal areal mass 20 g/cm2', 'equal20gcm2'))):
        en, T = wide_np(BB + f'wide_T_{kind}_{tag}.csv')
        rel = 100 * (T / np.nanmean(T, axis=1, keepdims=True) - 1)
        panel_id = 'ABCD'[row * 2 + col]
        _out(F, f'panel{panel_id}_{nm}_{lab.replace(" ", "_")}_relative_to_mean_percent', HDR,
             [[en[j]] + ['' if np.isnan(rel[j][i]) else '%.6g' % rel[j][i] for i in range(10)] for j in range(len(en))],
             f'panel {panel_id}: {ttl}, {lab}: transmission of each mixture relative to the mean over the ten mixtures (%)')
readme(F, f"""figs_new/fig_broadbeam_primary_vs_total.png - G-P (epixs_cli) exposure buildup, 2 x 2 panels, x log axis (Energy_MeV), y range -30..30, horizontal line at 0.
y = 100 x (T_mixture / mean over the ten mixtures - 1), % (already computed here, ready to plot). Empty cells = slab thinner than 0.5 mfp (equal-thickness 10 cm: E >= 3 MeV; 20 g/cm2: E >= 7 MeV), outside the G-P validity range 0.5-40 mfp, so no value exists (the PNG leaves them blank).
PANEL A (top-left):     equal thickness 10 cm, primary only          PANEL B (top-right):    equal areal mass 20 g/cm2, primary only
PANEL C (bottom-left):  equal thickness 10 cm, with buildup          PANEL D (bottom-right): equal areal mass 20 g/cm2, with buildup
The absolute transmissions these are computed from are in extra_no_figure/ (broadbeam_wide_T_*, broadbeam_wide_B_*).
CAVEAT: G-P accuracy for these low-Z mixtures is not independently verified (buildup cross-check was dropped); state the ranking claims as G-P-based.
{MIXNOTE}""")

F = 'fig_broadbeam_spread'
bs = list(csv.DictReader(open(BB + 'broadbeam_summary.csv')))
for pn, (kind, v, nm, ttl) in zip('AB', (('thickness_cm', '10', 'equal10cm', 'equal thickness 10 cm'), ('areal_mass_g_cm2', '20', 'equal20gcm2', 'equal areal mass 20 g/cm2'))):
    s = [x for x in bs if x['condition'] == kind and x['value'] == v]
    table(F, f'panel{pn}_{nm}_spread_percent', ['Energy_MeV', 'primary_only', 'with_buildup'], [[x['E_MeV'], x['spread_primary_pct'], x['spread_broad_pct']] for x in s],
          f'panel {pn}: {ttl}: best-to-worst mixture spread in transmission (%), primary only vs with buildup')
readme(F, """figs_new/fig_broadbeam_spread.png - G-P based, two panels, log-log axes, x = Energy_MeV.
PANEL A (left): equal thickness 10 cm.  PANEL B (right): equal areal mass 20 g/cm2.
y = best-to-worst mixture spread in transmission (%); two curves per panel: primary_only (circles) and with_buildup (squares).
The complete summary for all thicknesses / areal masses (1, 2, 5, 10, 20, 50, 100) is extra_no_figure/broadbeam_summary_all_conditions.csv.""")

# ---------------------------------------------------------------- data that is not (yet) in any figure + table
X = 'extra_no_figure'; os.makedirs(DEST + X, exist_ok=True)
for src, dst, what in ((SENS + 'element_mac.csv', 'element_MAC_pure_elements.csv', 'MAC of the pure elements C,H,N,O,S (XCOM) - basis of the sensitivity map, not plotted'),
                       (SENS + 'range_test_importance.csv', 'range_test_importance.csv', 'Pratt variance shares for the 7 range scenarios (robustness test of element importance), not plotted'),
                       (SENS + 'exchange_confirmation.csv', 'exchange_confirmation_full.csv', 'full exchange-confirmation table incl. additivity rows and absolute MACs'),
                       (SENS + 'summary.json', 'sensitivity_summary.json', 'summary numbers of the sensitivity analysis'),
                       (R + '/inputs/exchange_design.json', 'exchange_design.json', 'compositions used for the exchange runs'),
                       (BB + 'broadbeam_summary.csv', 'broadbeam_summary_all_conditions.csv', 'broad-beam spread/best/worst mixture for every thickness and areal mass'),
                       (BB + 'broadbeam_long.csv', 'broadbeam_long.csv', 'broad-beam per mixture, energy and condition (long format)')):
    shutil.copy(src, f'{DEST}{X}/{dst}'); INDEX.append((X, dst, what))
for fn in sorted(os.listdir(BB)):
    if fn.startswith('wide_'):
        shutil.copy(BB + fn, f'{DEST}{X}/broadbeam_{fn}'); INDEX.append((X, f'broadbeam_{fn}', 'absolute transmission or buildup factor (wide, energy x mixture); tags t=cm thickness, m=g/cm2 areal mass'))
open(f'{DEST}{X}/README.txt', 'w').write("""Data that no PNG in figs_new/ plots directly (kept for the paper text / possible new figures).
broadbeam_wide_T_primary_*/T_broad_*/B_* : absolute values behind fig_broadbeam_primary_vs_total (tags: t10cm/t20cm = equal thickness, m20gcm2/m50gcm2 = equal areal mass 20/50 g/cm2).
element_MAC_pure_elements.csv and range_test_importance.csv support the sensitivity figures and their robustness statement.
""")

# Table 2 data
table('.', 'Table2_mixture_composition_and_density', ['Mixture', 'C', 'H', 'O', 'N', 'S', 'Density_g_per_cm3', 'BSA_weight_fraction'],
      [[f'Mix{m}'] + ['%.6f' % COMP[m][e] for e in 'CHONS'] + ['%.6f' % RHO[m], '%.4f' % BSA[m]] for m in MIXES],
      'exact elemental mass fractions, density and BSA weight fraction of the ten mixtures (manuscript Table 2); not a figure')

# ---------------------------------------------------------------- index files
with open(DEST + 'README.txt', 'w') as f:
    f.write("""FIGURE DATA - recalculated CHON study (Geant4 + XCOM/G-P), regenerate with ~/CHON/recalc_2026-09-19/analysis/g4/export_figure_data.py
LAYOUT: one folder per figure. Folder name = PNG name in ../figs_new/ (e.g. fig_mac_lac/ <-> figs_new/fig_mac_lac.png).
Inside: panelA_..., panelB_... = panels in reading order (left->right, top->bottom); each folder has its own README.txt with the x column, y columns,
axis scales and what is a curve / error bar / shaded band. All files are CSV with header row; "wide" = Energy_MeV then one column per curve.
figure_data_all_sheets.xlsx = same tables, sheet 'INDEX' lists which sheet is which file. Mixture composition: Table2_mixture_composition_and_density.csv.
Not in any figure: fig1_geometry (schematic, no data).
Uncertainties are binomial from the pooled counts (2 x 10^6 photons per point); above 1 MeV the Geant4 composition spread is noise-limited (use XCOM there).

FIGURE FOLDER                              SOURCE                        CONTENT
""")
    src = {'fig_mac_lac': 'Geant4', 'fig_mfp_hvl_tvl': 'Geant4', 'fig_rpe': 'Geant4 (1 cm)', 'fig_validation_deviation': 'Geant4 vs XCOM', 'fig_seed_check': 'Geant4 (2 seeds)',
           'fig_zeff': 'XCOM', 'fig_zeff_baseline': 'XCOM', 'fig_process_fractions': 'XCOM (Mix1)', 'fig_ebf': 'G-P', 'fig_composition_sensitivity': 'Geant4 + XCOM',
           'fig_elemental_sensitivity': 'XCOM mixture rule', 'fig_library_vs_ensemble_correlation': 'XCOM ensemble', 'fig_element_importance': 'XCOM ensemble',
           'fig_exchange_confirmation': 'Geant4 vs XCOM', 'fig_broadbeam_primary_vs_total': 'G-P', 'fig_broadbeam_spread': 'G-P'}
    seen = []
    for fo, fi, what in INDEX:
        if fo not in seen:
            seen.append(fo)
            f.write(f"\n{fo + '/':<43}{src.get(fo, ''):<30}(see {fo}/README.txt)\n" if fo != '.' else '\n(top level, not a figure)\n')
        f.write(f"   {fi:<66}{what}\n")
from openpyxl import Workbook
wb = Workbook(); ws = wb.active; ws.title = 'INDEX'; ws.append(['sheet', 'figure folder', 'csv file', 'description'])
for fo, fi, what in INDEX:
    nm = _sheet(fo, fi[:-4]) if fo in ALIAS else None
    if nm in SHEETS: ws.append([nm, fo, fi, what])
for nm, rows in SHEETS.items():
    w = wb.create_sheet(nm); [w.append(r) for r in rows]
wb.save(DEST + 'figure_data_all_sheets.xlsx')
print(len(INDEX), 'files in', len({i[0] for i in INDEX}), 'folders written to', DEST)
