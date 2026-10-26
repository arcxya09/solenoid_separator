import tolB, sys, json, numpy as np
from tolB import *
def sim2(err=None,bk=None,dB=0.,seed=1,N=1500):
    wp.SRC=tolB.SRC; nm='tb_%g'%dB; mk30(nm,tolB.B0*(1+dB)); sol=get_sol(nm); Bf=ErrField(sol,r0=0.104,**err) if err else sol
    b=dict(dR=0.0); b.update(bk or {}); zp=np.arange(tolB.ZD-0.06,tolB.ZD+0.06,0.005)
    recs,zp,_=run_beam(nm,tolB.SRC,N=N,B=Bf,I=tolB.I*1e-3,zplanes=zp,seed=seed,beam_kw=b,iters=5); return recs,zp
tolB.sim=sim2
def run(case):
    I,src=case; tolB.I=I; tolB.SRC=src
    with Pool(4) as p: mc=p.map(mcjob,range(20),chunksize=1)
    for lev in ('L0','L1'):
        T=np.array([r[lev]['T'] for r in mc]); P=np.array([r[lev]['pur'] for r in mc])
        print('I=%g'%I,lev,'mean T %.4f P10 %.4f min %.4f minpur %.5f yield %.2f'%(T.mean(),np.percentile(T,10),T.min(),P.min(),((T>=0.99)&(P>=0.999)).mean()),flush=True)
    return mc
if __name__=='__main__':
    out={}
    for case in ((6.0,'srcB_6mA_V17400.pkl'),(5.0,'src104_30k_5mA.pkl')): out[case[0]]=run(case)
    json.dump(dict(note='b2,a2 defined at r0=0.104 m (beam radius)',SIG=tolB.SIG,mc=out),open('tolB2.json','w'))
