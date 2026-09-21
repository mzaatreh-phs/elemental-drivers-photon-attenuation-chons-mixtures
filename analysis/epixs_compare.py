#!/usr/bin/env python3
"""Compare epixs_cli G-P (mixture/Zeq path) against the EpiXS exports for the ten ORIGINAL CHON mixtures.
Exports: ~/CHON/EXPIXS DATA/{Exposure,Energy Absorption} Buildup Factor -MIXnn.csv   (Standard grid rows)."""
import sys,os,csv,glob,statistics as st
sys.path.insert(0,os.path.expanduser('~/epixs_cli')); sys.path.insert(0,os.path.expanduser('~/epixs_cli/epixs'))
from zeff import zeq
from gp_calculator import material_gp_coeffs, gp_buildup_factor
OLD={1:(.5166,.0745,.3281,.0808),2:(.5531,.0795,.2867,.0808),3:(.4802,.0695,.3696,.0808),4:(.5013,.0675,.3302,.1010),5:(.5439,.0771,.2498,.1293),6:(.5257,.0754,.3020,.0970),7:(.5075,.0737,.3542,.0646),8:(.4984,.0728,.3803,.0485),9:(.4893,.0720,.4064,.0323),10:(.4712,.0703,.4587,0.0)}
MFP=[1,2,3,4,5,6,7,8,9,10,11,12,15,20,25,30,35,40]
D=os.path.expanduser('~/CHON/archive/CHON_root/EXPIXS DATA/')
out=[];summary={}
for kind,label in (("exposure","Exposure Buildup Factor -MIX%02d.csv"),("absorption","Energy Absorption Buildup Factor -MIX%02d.csv")):
    for m,(c,h,o,n) in OLD.items():
        comp={k:v for k,v in dict(C=c,H=h,O=o,N=n).items() if v>0}
        rows=[r for r in csv.reader(open(D+label%m,encoding='utf-8-sig')) if r and r[0]=='Standard']
        for r in rows:
            E=float(r[1])/1000.0; zeq_ex=float(r[3]); vals=[float(x) for x in r[9:9+18]]
            try:
                z,z1,z2=zeq(comp,E)
                b,cc,a,Xk,d=material_gp_coeffs(E,z,z1,z2,kind)
                mine=[gp_buildup_factor(x,b,cc,a,Xk,d) for x in MFP]
            except Exception as ex:
                out.append((kind,m,E,zeq_ex,None,None,None,str(ex)[:60])); continue
            rel=[100*(p/q-1) for p,q in zip(mine,vals)]
            out.append((kind,m,E,zeq_ex,z,max(abs(x) for x in rel),rel[MFP.index(40)],""))
with open(os.path.expanduser('~/CHON/recalc_2026-09-19/analysis/epixs_compare.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(["kind","mix","E_MeV","Zeq_EpiXS","Zeq_epixs_cli","max_abs_rel_diff_%_1to40mfp","rel_diff_%_at_40mfp","error"]); w.writerows(out)
for kind in ("exposure","absorption"):
    ok=[r for r in out if r[0]==kind and r[4] is not None]; bad=[r for r in out if r[0]==kind and r[4] is None]
    mx=[r[5] for r in ok]; dz=[abs(r[4]/r[3]-1)*100 for r in ok]
    print(f"\n== {kind}: {len(ok)} (mix,E) rows compared, {len(bad)} failed")
    print(f"   max |diff| over all depths: median {st.median(mx):.3f}%  95th pct {sorted(mx)[int(.95*len(mx))]:.3f}%  max {max(mx):.3f}%")
    print(f"   Zeq relative difference:   median {st.median(dz):.3f}%  max {max(dz):.3f}%")
    print("   by energy (max over mixtures of max|diff| %):")
    for E in sorted({r[2] for r in ok}):
        v=[r[5] for r in ok if r[2]==E]; print(f"     {E:8.3f} MeV  max {max(v):8.3f}%  median {st.median(v):8.3f}%")
    if bad: print("   failures:",bad[:3])
