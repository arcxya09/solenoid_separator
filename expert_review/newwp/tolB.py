from hc30 import *
import hc30
SIG=dict(coil_dx=0.5e-3,coil_tilt=0.5e-3,src_dx=0.5e-3,src_tilt=1.0e-3,src_dR=0.25e-3,dB=5e-4,dE=5e-4,b2=2.5e-4,col_dx=0.3e-3,col_dz=5e-3)
B0=1.912; I=6.0; SRC='srcB_6mA_V17400.pkl'; ZD=2.415; A=0.021
def sim(err=None,bk=None,dB=0.,seed=1,N=1500):
    wp.SRC=SRC; nm='tb_%g'%dB; mk30(nm,B0*(1+dB)); sol=get_sol(nm); Bf=ErrField(sol,r0=0.03,**err) if err else sol
    b=dict(dR=0.0); b.update(bk or {})
    zp=np.arange(ZD-0.06,ZD+0.06,0.005)
    recs,zp,_=run_beam(nm,SRC,N=N,B=Bf,I=I*1e-3,zplanes=zp,seed=seed,beam_kw=b,iters=5)
    return recs,zp
def ev(recs,zp,cx=0.,cy=0.,dz=0.):
    L0=fixed(recs,zp,ZD+dz,A,cx=cx,cy=cy); w=abs(zp-ZD)<=0.05; L1=stage({k:v[w] for k,v in recs.items()},zp[w],A); return L0,L1
def oat(a):
    par,amp=a; s=SIG[par]*amp; err={}; bk={}; dB=0.
    if par=='coil_dx': err=dict(dx=s)
    elif par=='coil_tilt': err=dict(tx=s)
    elif par=='b2': err=dict(b2=s)
    elif par=='src_dx': bk=dict(dx=s)
    elif par=='src_tilt': bk=dict(tx=s)
    elif par=='src_dR': bk=dict(dR=s)
    elif par=='dE': bk=dict(dE=s)
    elif par=='dB': dB=s
    L0,L1=ev(*sim(err,bk,dB)); return dict(par=par,amp=amp,L0=L0,L1=L1)
def mcjob(i):
    rng=np.random.default_rng(2000+i); g=lambda key: SIG[key]*rng.standard_normal()
    err=dict(dx=g('coil_dx'),dy=g('coil_dx'),tx=g('coil_tilt'),ty=g('coil_tilt'),b2=g('b2'),a2=g('b2'))
    bk=dict(dx=g('src_dx'),dy=g('src_dx'),tx=g('src_tilt'),ty=g('src_tilt'),dR=g('src_dR'),dE=g('dE'))
    recs,zp=sim(err,bk,g('dB'),seed=i); L0,L1=ev(recs,zp,g('col_dx'),g('col_dx'),g('col_dz')); return dict(i=i,L0=L0,L1=L1)
if __name__=='__main__':
    with Pool(4) as p:
        nom=p.map(oat,[('dB',0)])[0]; print('nominal',nom,flush=True)
        J=[(par,amp) for par in ('dB','dE','src_tilt','src_dx','src_dR','coil_tilt','coil_dx','b2') for amp in (2,4)]
        oa=p.map(oat,J,chunksize=1)
        for r in oa: print('%-9s x%d  L0 T=%.3f pur=%.5f | L1 T=%.3f pur=%.5f z=%.3f'%(r['par'],r['amp'],r['L0']['T'],r['L0']['pur'],r['L1']['T'],r['L1']['pur'],r['L1']['z']),flush=True)
        mc=p.map(mcjob,range(20),chunksize=1)
    json.dump(dict(SIG=SIG,B0=B0,I=I,ZD=ZD,A=A,nom=nom,oat=oa,mc=mc),open('tolB.json','w'))
    for lev in ('L0','L1'):
        T=np.array([r[lev]['T'] for r in mc]); P=np.array([r[lev]['pur'] for r in mc])
        print(lev,'mean T %.4f P10 %.4f min %.4f minpur %.5f yield %.2f'%(T.mean(),np.percentile(T,10),T.min(),P.min(),((T>=0.99)&(P>=0.999)).mean()),flush=True)
