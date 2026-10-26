# one-at-a-time sensitivities + multi-error MC for the old and new working points
import wp
from wp import *
from multiprocessing import Pool
SIG=dict(coil_dx=0.5e-3,coil_tilt=0.5e-3,src_dx=0.5e-3,src_tilt=1.5e-3,src_dR=0.25e-3,dB=5e-4,dE=5e-4,b2=5e-4,col_dx=0.3e-3,col_dz=5e-3)
CASES={  # name: R0, B, I(mA), src
 'new0':(0.06,1.894,0.0,'src60_0p2.pkl'), 'new05':(0.06,1.9129,0.5,'src60_0p5.pkl'),
 'base0':(0.03,2.800,0.0,'src_v3.pkl'),  'base015':(0.03,2.838,0.15,'src_v3.pkl')}
def setup(case,dB=0.):
    R0,B,I,src=CASES[case]; wp.SRC=src; nm='%s_%g'%(case,dB); mk(nm,R0,B*(1+dB),'C3'); return nm,I
def sim(case,err=None,bk=None,dB=0.,N=1500,seed=1):
    nm,I=setup(case,dB); d=DESIGNS[nm]; sol=get_sol(nm); Bf=ErrField(sol,**err) if err else sol
    Rs=round(float(np.median(load_src(wp.SRC)['r'])),3); b=dict(dR=d['R']-Rs); b.update(bk or {})
    if 'dR' in (bk or {}): b['dR']=d['R']-Rs+bk['dR']
    zp=np.arange(1.30,1.43,0.0025)
    recs,zp,_=run_beam(nm,wp.SRC,N=N,B=Bf,I=I*1e-3,zplanes=zp,seed=seed,beam_kw=b,iters=4)
    return recs,zp
def nominal(case):
    recs,zp=sim(case); k=int(np.argmin(abs(zp-1.365))); return at_plane(recs,zp,k,0.999,recenter=False)
def oat(a):
    case,par,amp,A=a; s=SIG[par]*amp; err={}; bk={}; dB=0.
    if par=='coil_dx': err=dict(dx=s)
    elif par=='coil_tilt': err=dict(tx=s)
    elif par=='b2': err=dict(b2=s)
    elif par=='src_dx': bk=dict(dx=s)
    elif par=='src_tilt': bk=dict(tx=s)
    elif par=='src_dR': bk=dict(dR=s)
    elif par=='dE': bk=dict(dE=s)
    elif par=='dB': dB=s
    recs,zp=sim(case,err,bk,dB)
    L0=fixed(recs,zp,1.365,A); w=abs(zp-1.365)<=0.05; L1=stage({k:v[w] for k,v in recs.items()},zp[w],A)
    return dict(case=case,par=par,amp=amp,L0=L0,L1=L1)
def mcjob(a):
    case,i,k,A=a; rng=np.random.default_rng(1000+i); g=lambda key: SIG[key]*k*rng.standard_normal()
    err=dict(dx=g('coil_dx'),dy=g('coil_dx'),tx=g('coil_tilt'),ty=g('coil_tilt'),b2=g('b2'),a2=g('b2'))
    bk=dict(dx=g('src_dx'),dy=g('src_dx'),tx=g('src_tilt'),ty=g('src_tilt'),dR=g('src_dR'),dE=g('dE'))
    recs,zp=sim(case,err,bk,g('dB'),seed=i)
    cx,cy,dz=g('col_dx'),g('col_dx'),g('col_dz')
    L0=fixed(recs,zp,1.365+dz,A,cx=cx,cy=cy); w=abs(zp-1.365)<=0.05; L1=stage({kk:v[w] for kk,v in recs.items()},zp[w],A)
    return dict(case=case,i=i,k=k,L0=L0,L1=L1)
if __name__=='__main__':
    with Pool(4) as p: nom=dict(zip(CASES,p.map(nominal,list(CASES))))
    for c,v in nom.items(): print('nominal',c,v,flush=True)
    A={c:v['a'] for c,v in nom.items()}
    J=[(c,par,amp,A[c]) for c in CASES for par in ('dB','dE','src_tilt','src_dx','src_dR','coil_tilt','coil_dx','b2') for amp in (2,4)]
    with Pool(4) as p: oa=p.map(oat,J,chunksize=1)
    for r in oa: print('%-8s %-9s x%d  L0 T=%.3f pur=%.5f | L1 T=%.3f'%(r['case'],r['par'],r['amp'],r['L0']['T'],r['L0']['pur'],r['L1']['T']),flush=True)
    J=[(c,i,1.0,A[c]) for c in CASES for i in range(20)]
    with Pool(4) as p: mc=p.map(mcjob,J,chunksize=1)
    json.dump(dict(SIG=SIG,nom=nom,oat=oa,mc=mc),open('tol.json','w'))
    for c in CASES:
        for lev in ('L0','L1'):
            T=np.array([r[lev]['T'] for r in mc if r['case']==c]); P=np.array([r[lev].get('pur',1.0) for r in mc if r['case']==c])
            print(c,lev,'mean T %.4f P10 %.4f min %.4f minpur %.5f yield %.2f'%(T.mean(),np.percentile(T,10),T.min(),P.min(),((T>=0.99)&(P>=0.999)).mean()),flush=True)
