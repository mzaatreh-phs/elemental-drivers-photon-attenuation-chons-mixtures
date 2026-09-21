#!/usr/bin/env python3
"""XCOM/G-P quantities for the exact CHONS compositions (epixs_cli): Zeff (direct), process fractions + crossover, EBF/EABF."""
import sys,os,csv,json
sys.path.insert(0,os.path.expanduser('~/epixs_cli')); sys.path.insert(0,os.path.expanduser('~/epixs_cli/epixs'))
from zeff import zeff, zeq
from gp_calculator import material_gp_coeffs, gp_buildup_factor
from epixs.xcom_engine import mass_attenuation
R=os.path.expanduser('~/CHON/recalc_2026-09-19'); O=R+'/analysis/xcom/'
inp=json.load(open(R+'/inputs/inputs_exact.json')); mix=inp['mixtures']
E38=sorted(inp['production_energies']+inp['validation_energies'])
MFP=[1,2,3,4,5,6,7,8,9,10,11,12,15,20,25,30,35,40]
OLD={1:(.5166,.0745,.3281,.0808),2:(.5531,.0795,.2867,.0808),3:(.4802,.0695,.3696,.0808),4:(.5013,.0675,.3302,.1010),5:(.5439,.0771,.2498,.1293),6:(.5257,.0754,.3020,.0970),7:(.5075,.0737,.3542,.0646),8:(.4984,.0728,.3803,.0485),9:(.4893,.0720,.4064,.0323),10:(.4712,.0703,.4587,0.0)}
comps={('new',int(m)):x['composition'] for m,x in mix.items()}
comps.update({('old',m):{k:v for k,v in dict(C=a,H=b,O=c,N=d).items() if v>0} for m,(a,b,c,d) in OLD.items()})
refs={'water':{'H':0.111894,'O':0.888106},'PMMA':{'C':0.599848,'H':0.080538,'O':0.319614}}
# 1) Zeff
with open(O+'zeff.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['set','material','E_MeV','Zeff_direct'])
    for (s,m),c in comps.items():
        for E in E38: w.writerow([s,f'Mix{m}',E,f'{zeff(c,E):.6f}'])
    for n,c in refs.items():
        for E in E38: w.writerow(['ref',n,E,f'{zeff(c,E):.6f}'])
print('zeff done')
# 2) process fractions + crossover (XCOM components)
fine=[10**(1+1.0*i/399) /1000.0 for i in range(400)]   # 10..100 keV, log spaced (MeV)
with open(O+'process_fractions.csv','w',newline='') as f, open(O+'crossover.csv','w',newline='') as g:
    w=csv.writer(f); w.writerow(['set','material','E_MeV','F_photoelectric','F_incoherent','F_coherent','F_pair','mac_total_cm2_g'])
    x=csv.writer(g); x.writerow(['set','material','crossover_keV'])
    for (s,m),c in comps.items():
        r=mass_attenuation(c,E38); tot=r['total_with_coherent']
        for i,E in enumerate(E38):
            pair=r['pair_nuclear'][i]+r['pair_electron'][i]
            w.writerow([s,f'Mix{m}',E]+[f'{v/tot[i]:.6f}' for v in (r['photoelectric'][i],r['incoherent'][i],r['coherent'][i],pair)]+[f'{tot[i]:.6f}'])
        rf=mass_attenuation(c,fine); d=[p-q for p,q in zip(rf['photoelectric'],rf['incoherent'])]
        cross=None
        for i in range(len(d)-1):
            if d[i]>0>=d[i+1]:
                import math
                t=d[i]/(d[i]-d[i+1]); cross=1000*math.exp(math.log(fine[i])+t*(math.log(fine[i+1])-math.log(fine[i]))); break
        x.writerow([s,f'Mix{m}',f'{cross:.3f}' if cross else 'none'])
print('process fractions done')
# 3) EBF / EABF (G-P via Zeq), E >= 0.015 MeV
E37=[E for E in E38 if E>=0.015]
for kind in ('exposure','absorption'):
    with open(O+f'{kind}_buildup.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['set','material','E_MeV','Zeq']+[f'B_{x}mfp' for x in MFP])
        for (s,m),c in comps.items():
            for E in E37:
                z,z1,z2=zeq(c,E); b,cc,a,Xk,d=material_gp_coeffs(E,z,z1,z2,kind)
                w.writerow([s,f'Mix{m}',E,f'{z:.5f}']+[f'{gp_buildup_factor(x,b,cc,a,Xk,d):.6g}' for x in MFP])
    print(kind,'done')
