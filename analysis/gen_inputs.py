#!/usr/bin/env python3
"""Generate exact inputs and Geant4 (brks) macros for the CHON recalculation.
Compositions: exact formulas (BSA from UniProt P02769 incl. S).  Densities: additive-volume rule.
Thickness per energy: optical depth 1.5 from XCOM (epixs_cli) for the exact composition."""
import json,os,sys
sys.path.insert(0,os.path.expanduser('~/epixs_cli'))
from epixs.xcom_engine import mass_attenuation
R=os.path.expanduser('~/CHON/recalc_2026-09-19')
d=json.load(open('/tmp/claude-1000/-home-mzaatreh/ca0288cb-0500-48f4-9790-a869643806ea/scratchpad/comp.json'))
RHO=dict(cellulose=1.50,BSA=1.33,sorbitol=1.49,stearic=0.94)          # g/cm3, as cited in the manuscript
T1={1:(30,50,10,10),2:(30,50,0,20),3:(30,50,20,0),4:(37.5,62.5,0,0),5:(0,80,10,10),6:(20,60,10,10),7:(40,40,10,10),8:(50,30,10,10),9:(60,20,10,10),10:(80,0,10,10)}
PROD=[0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
VAL=[0.356,0.511,0.662,1.173,1.33,2.51]
GRID=sorted(PROD+VAL)
OD=1.5
mix={}
for m,ex,old in d["rows"]:
    w=dict(zip(("cellulose","BSA","sorbitol","stearic"),[x/100 for x in T1[m]]))
    rho=1/sum(w[k]/RHO[k] for k in w)
    comp={e:round(ex[e],10) for e in "CHONS" if ex[e]>0}
    mac=mass_attenuation(comp,GRID)['total_with_coherent']
    thick=[OD/(x*rho) for x in mac]
    mix[m]=dict(composition=comp,density=rho,weights=w,energies=GRID,xcom_mac=mac,thickness_cm=thick)
json.dump(dict(rho_constituents=RHO,optical_depth_target=OD,production_energies=PROD,validation_energies=VAL,mixtures=mix),open(f'{R}/inputs/inputs_exact.json','w'),indent=1)
def macro(m,energies,seed,events,tag):
    x=mix[m]; idx=[GRID.index(e) for e in energies]
    L=["# CHON recalculation %s  Mix_%d seed=%d events=%d"%(tag,m,seed,events),
       "/control/verbose 1","/run/verbose 1",
       "/random/setSeeds %d %d"%(seed,seed+1),
       "/det/clearMaterial","/det/setMaterialName Mix_%d"%m,"/det/setDensity %.6f g/cm3"%x['density'],
       "/det/setThickness 1.0 cm"]
    L+=["/det/addElement %s %.10f"%(e,f) for e,f in x['composition'].items()]
    L+=["/run/setCut 0.001 mm","/det/createMaterial",
        "/scan/setEnergyGrid "+",".join(str(e) for e in energies),
        "/scan/setThicknessGrid "+",".join("%.6f"%x['thickness_cm'][i] for i in idx),
        "/scan/setEvents %d"%events,"/run/initialize"]
    return "\n".join(L)+"\n"
os.makedirs(f'{R}/macros/pilot',exist_ok=True); os.makedirs(f'{R}/macros/full',exist_ok=True)
PE=[0.01,0.03,0.1,0.662,1.173,2.51,10,15]
for m in (1,5,10):
    open(f'{R}/macros/pilot/Mix_{m}_pilot.mac','w').write(macro(m,PE,20260919,1000000,'PILOT'))
for m in mix:
    for k,seed in enumerate((1910001,1910101),1):
        open(f'{R}/macros/full/Mix_{m}_s{k}.mac','w').write(macro(m,GRID,seed,1000000,'FULL'))
for m in (1,5,10): print(m,"rho=%.5f"%mix[m]['density'],{e:round(f,4) for e,f in mix[m]['composition'].items()})
print("thickness Mix_1 (cm) 0.01/0.1/1/15 MeV:",[round(mix[1]['thickness_cm'][GRID.index(e)],2) for e in (0.01,0.1,1,15)])
