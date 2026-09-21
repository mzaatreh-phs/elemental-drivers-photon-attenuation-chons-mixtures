#!/usr/bin/env python3
import os, csv, math, numpy as np
A = os.path.expanduser('~/CHON/recalc_2026-09-19/analysis/')
def rows(p): return list(csv.DictReader(open(A + p)))
EB = [r for r in rows('xcom/exposure_buildup.csv') if r['set'] == 'new']
Es = sorted({float(r['E_MeV']) for r in EB}); print('EBF energies', len(Es), Es[0], Es[-1]); print('cols', [k for k in EB[0].keys()][:8])
def B(m, e, d): return float([r for r in EB if r['material'] == 'Mix%d' % m and abs(float(r['E_MeV']) - e) < 1e-9][0]['B_%dmfp' % d])
M = range(1, 11)
for d in (5, 10, 40):
    pk = [(max(((B(m, e, d), e) for e in Es))) for m in M]; print('depth %2d mfp: peak EBF range %.4g-%.4g at E=%s' % (d, min(p[0] for p in pk), max(p[0] for p in pk), sorted({p[1] for p in pk})))
print('== 40 mfp EBF ranges at low energies and ranking')
for e in (0.015, 0.02, 0.03, 0.04, 0.05, 0.08, 0.1, 0.15, 0.3, 1, 5, 10, 15):
    v = np.array([B(m, e, 40) for m in M]); print('  E=%-5g EBF(40) %.4g (Mix%d) - %.4g (Mix%d); spread %.1f%%; at 5 mfp spread %.1f%%' % (e, v.min(), v.argmin() + 1, v.max(), v.argmax() + 1, 100 * (v.max() / v.min() - 1), 100 * (lambda w: w.max() / w.min() - 1)(np.array([B(m, e, 5) for m in M]))))
print('  40 mfp at 0.1 MeV: Mix5 %.0f Mix10 %.0f Mix2 %.0f Mix4 %.0f Mix3 %.0f Mix6 %.0f' % tuple(B(m, 0.1, 40) for m in (5, 10, 2, 4, 3, 6)))
print('  total-response proxy exp(-40)*B at 0.1 MeV: min %.3g (Mix%d) max %.3g (Mix%d)' % (math.exp(-40) * min(B(m, .1, 40) for m in M), min(M, key=lambda m: B(m, .1, 40)), math.exp(-40) * max(B(m, .1, 40) for m in M), max(M, key=lambda m: B(m, .1, 40))))
print('== rank reversals: highest/lowest EBF mixture at 40 mfp per energy (changes only)')
last = None
for e in Es:
    v = np.array([B(m, e, 40) for m in M]); o = (v.argmax() + 1, v.argmin() + 1)
    if o != last: print('   from E=%g: highest Mix%d lowest Mix%d' % (e, *o)); last = o
print('== broad-beam summary')
S = rows('broadbeam/broadbeam_summary.csv'); print(list(S[0].keys()))
for kind, val in (('thickness_cm', '10'), ('areal_mass_g_cm2', '20')):
    s = [r for r in S if r['condition'] == kind and r['value'] == val]
    print(' %s=%s: %d energies (%g-%g MeV)' % (kind, val, len(s), float(s[0]['E_MeV']), float(s[-1]['E_MeV'])))
    for r in s:
        if float(r['E_MeV']) in (0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1, 1.5, 2, 3): print('   E=%-5s spread primary %.2f%% total %.2f%% | best(primary) Mix%s best(total) Mix%s worst(total) Mix%s | Spearman %s' % (r['E_MeV'], float(r['spread_primary_pct']), float(r['spread_broad_pct']), r['best_mix_primary'], r['best_mix_broad'], r['worst_mix_broad'], r['spearman_primary_vs_broad'][:5]))
    diff = sum(r['best_mix_primary'] != r['best_mix_broad'] for r in s); print('   best mixture differs (primary vs total) at %d of %d energies' % (diff, len(s)))
    comp = [(float(r['spread_primary_pct']) / float(r['spread_broad_pct']), float(r['E_MeV'])) for r in s if 0.1 <= float(r['E_MeV']) <= 1.0]; print('   spread ratio primary/total for 0.1-1 MeV: %.1f - %.1f' % (min(c[0] for c in comp), max(c[0] for c in comp)))
    from itertools import groupby
    bp = [(float(r['E_MeV']), r['best_mix_primary'], r['best_mix_broad']) for r in s]; print('   best mix by energy:', bp[:8], '...')
    sp = [float(r['spearman_primary_vs_broad']) for r in s if r['spearman_primary_vs_broad'] not in ('', 'nan')]; print('   Spearman range %.2f..%.2f; median %.2f' % (min(sp), max(sp), np.median(sp)))
    hi = [(float(r['E_MeV']), float(r['spread_primary_pct']), float(r['spread_broad_pct'])) for r in s if float(r['E_MeV']) >= 1.3]; print('   E>=1.3 MeV spreads max primary %.2f total %.2f' % (max(h[1] for h in hi), max(h[2] for h in hi)))
