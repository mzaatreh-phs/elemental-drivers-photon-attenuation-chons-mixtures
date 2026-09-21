#!/usr/bin/env python3
"""Designed controlled-exchange set for Geant4 confirmation of the mixture-rule sensitivities.
Base = mean library composition and density; each variant replaces carbon by one element (fixed density, fixed 1.0 cm slab)."""
import json, os
import numpy as np
R = os.path.expanduser('~/CHON/recalc_2026-09-19')
inp = json.load(open(R + '/inputs/inputs_exact.json')); mix = inp['mixtures']
EL = ['C', 'H', 'N', 'O', 'S']
W = np.array([[mix[str(m)]['composition'].get(e, 0.0) for e in EL] for m in range(1, 11)])
base = W.mean(0); rho = float(np.mean([mix[str(m)]['density'] for m in range(1, 11)]))
def variant(d):  # d = {element: +delta}, carbon takes the balance
    w = base.copy()
    for e, x in d.items(): w[EL.index(e)] += x; w[0] -= x
    return w
DES = {'base': {}, 'dO': {'O': 0.05}, 'dN': {'N': 0.05}, 'dS': {'S': 0.015}, 'dH': {'H': 0.015}, 'dALL': {'O': 0.05, 'N': 0.05, 'S': 0.015, 'H': 0.015}}
ENERGIES = [0.01, 0.02, 0.03, 0.05, 0.1, 0.3, 1.0, 10.0]
design = {k: {'composition': dict(zip(EL, variant(v).round(10).tolist())), 'delta_wt': v} for k, v in DES.items()}
json.dump({'density_g_cm3': rho, 'energies_MeV': ENERGIES, 'design': design}, open(R + '/inputs/exchange_design.json', 'w'), indent=1)
jobs = []
for i, E in enumerate(ENERGIES):
    for si, seed_base in enumerate((3000000, 4000000), 1):
        for k, ci in zip(design, range(len(design))):
            name = f's{si}_{k}_E{i:02d}'; seed = seed_base + 1000 * ci + 10 * i + 7
            L = [f'# CHON controlled exchange {k} E={E} MeV', '/control/verbose 1', '/run/verbose 1', f'/random/setSeeds {seed} {seed + 1}', '/det/clearMaterial',
                 f'/det/setMaterialName {k}', f'/det/setDensity {rho:.6f} g/cm3', '/det/setThickness 1.0 cm']
            L += [f'/det/addElement {e} {f:.10f}' for e, f in design[k]['composition'].items() if f > 0]
            L += ['/run/setCut 0.001 mm', '/det/createMaterial', f'/scan/setEnergyGrid {E}', '/scan/setEvents 1000000', '/run/initialize']
            open(f'{R}/macros/exchange/{name}.mac', 'w').write('\n'.join(L) + '\n'); jobs.append(name)
open(R + '/macros/exchange/jobs.txt', 'w').write('\n'.join(jobs) + '\n')
print(len(jobs), 'jobs; base composition', dict(zip(EL, base.round(5))), 'rho %.4f' % rho)
for k, v in design.items(): print(f'  {k:5s}', {e: round(f, 5) for e, f in v['composition'].items()}, 'sum=%.6f' % sum(v['composition'].values()))
