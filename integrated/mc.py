# Multi-error Monte Carlo: python3 mc.py design src sigma_scale nsamp out.json [I f0]
# every sample draws all errors at once; evaluated (L0) with the collector fixed at its nominal position
# (plus its own placement error) and (L1) with the collector x/y/z stage re-optimised on the Rb-87 spot.
from sep import *
import sys
SIG=dict(coil_dx=0.5e-3,coil_tilt=0.5e-3,src_dx=0.5e-3,src_tilt=1.5e-3,src_dR=0.25e-3,
         dB=3e-4,dE=3e-4,b2=5e-4,col_dx=0.3e-3,col_dz=5e-3,einzel=100.0)
if os.environ.get('SIGJSON'): SIG.update(json.loads(os.environ['SIGJSON']))
def draw(rng,k):
    g=lambda key: SIG[key]*k*rng.standard_normal()
    return dict(err=dict(dx=g('coil_dx'),dy=g('coil_dx'),tx=g('coil_tilt'),ty=g('coil_tilt'),b2=g('b2'),a2=g('b2')),
                beam_kw=dict(dx=g('src_dx'),dy=g('src_dx'),tx=g('src_tilt'),ty=g('src_tilt'),dR=g('src_dR'),dE=g('dE')),
                dB=g('dB'),col=(g('col_dx'),g('col_dx'),g('col_dz')),dVl=g('einzel'))
def job(a):
    design,src,src2,k,i,N,I,f0=a; d=DESIGNS[design]; rng=np.random.default_rng(1000+i)
    e=draw(rng,k); sol=get_sol(design,d['B0']*(1+e['dB'])); B=ErrField(sol,**e['err'],r0=d['R'])   # b2/a2 referred to the beam radius R0
    bk=dict(e['beam_kw'])
    if src2: bk.update(src2=src2,wmix=e['dVl']/float(os.environ.get('SRC2DV','1000')))       # src2 = source run at Vl-1 kV (more negative)
    zp=np.arange(d['zc']-0.06,d['zc']+0.06,0.0025)
    recs,zp,_=run_beam(design,src,N=N,B=B,I=I,fneu=fneu_profile(design,f0) if f0 else None,zplanes=zp,seed=i,beam_kw=bk,iters=int(os.environ.get('ITERS','4')))
    cx,cy,dz=e['col']
    L0=fixed(recs,zp,d['zc']+dz,d['a'],cx=cx,cy=cy)
    w=abs(zp-d['zc'])<=0.05; L1=stage({kk:v[w] for kk,v in recs.items()},zp[w],d['a'])
    bb=best(recs,zp)
    e['err']={kk:float(v) for kk,v in e['err'].items()}; e['beam_kw']={kk:float(v) for kk,v in e['beam_kw'].items()}
    return dict(i=i,e=e,L0=L0,L1=L1,best=bb)
if __name__=='__main__':
    from multiprocessing import Pool
    design,src=sys.argv[1],sys.argv[2]; k=float(sys.argv[3]); n=int(sys.argv[4]); out=sys.argv[5]
    I=float(sys.argv[6]) if len(sys.argv)>6 else 0.0; f0=float(sys.argv[7]) if len(sys.argv)>7 else 0.0
    src2=os.environ.get('SRC2'); N=int(os.environ.get('NMC','3000'))
    with Pool(4) as p: res=p.map(job,[(design,src,src2,k,i,N,I,f0) for i in range(n)],chunksize=1)
    json.dump(dict(SIG=SIG,k=k,I=I,f0=f0,res=res),open(out,'w'))
    for lev in ('L0','L1'):
        T=np.array([r[lev]['T'] for r in res]); P=np.array([r[lev].get('pur',1.0) for r in res])
        ok=(T>=0.99)&(P>=0.999)
        print(design,'k=%g'%k,lev,'mean T %.4f  P10 T %.4f  min T %.4f  min pur %.5f  yield(T>=99%%,pur>=99.9%%) %.3f'%(T.mean(),np.percentile(T,10),T.min(),P.min(),ok.mean()))
