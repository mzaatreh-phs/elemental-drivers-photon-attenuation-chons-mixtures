#!/usr/bin/env python3
"""One Geant4 job per (density label, mixture, energy): fixed 1.0 cm slab, 1e6 photons, one seed each.
Usage: gen_jobs.py densities.json    where densities.json = {"label": {"1": rho1, ..., "10": rho10}, ...}"""
import json,os,sys
R=os.path.expanduser('~/CHON/recalc_2026-09-19')
inp=json.load(open(f'{R}/inputs/inputs_exact.json'))
mix=inp['mixtures']; GRID=sorted(inp['production_energies']+inp['validation_energies'])
dens=json.load(open(sys.argv[1])); labels=list(dens)
os.makedirs(f'{R}/macros/points',exist_ok=True)
jobs=[]
for idx,E in enumerate(GRID):                      # outer loop = energy so all mixtures advance together
    for li,lab in enumerate(labels):
        for m in mix:
            name=f'{lab}_Mix_{m}_E{idx:02d}'
            seed=1_000_000*(li+1)+10_000*int(m)+idx*100+7
            L=[f'# CHON recalculation  label={lab} Mix_{m} E={E} MeV  1.0 cm fixed  events=1e6',
               '/control/verbose 1','/run/verbose 1',f'/random/setSeeds {seed} {seed+1}',
               '/det/clearMaterial',f'/det/setMaterialName Mix_{m}',f'/det/setDensity {dens[lab][m]:.6f} g/cm3','/det/setThickness 1.0 cm']
            L+=[f'/det/addElement {e} {f:.10f}' for e,f in mix[m]['composition'].items()]
            L+=['/run/setCut 0.001 mm','/det/createMaterial',f'/scan/setEnergyGrid {E}','/scan/setEvents 1000000','/run/initialize']
            open(f'{R}/macros/points/{name}.mac','w').write('\n'.join(L)+'\n'); jobs.append(name)
open(f'{R}/macros/points/jobs.txt','w').write('\n'.join(jobs)+'\n')
print(len(jobs),'jobs;',len(labels),'density sets;',len(GRID),'energies;',len(mix),'mixtures')
