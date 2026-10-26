# scaled design: E=30 keV, R=50 mm, coil L=1.7 m, bore 150-180 mm
from ring import *
import ring, fromsrc, pickle, json, sys, os
from scan import metrics
KW=dict(L=1.7,a_in=0.15,a_out=0.18,nz=340,nr=3,rmax=0.13,zmin=-1.0,zmax=5.0,dr=1e-3,dz=2.5e-3)
E=30e3; R=0.05; zs=float(os.environ.get('ZS','-0.6'))
def sol(B0): return get_sol(B0,**KW)
def beam(s,M,N,rng,scale_ang=None):
    if scale_ang is None: scale_ang=float(os.environ.get('SCALEANG',str(1/np.sqrt(3))))
    r,vr,vt,vz=fromsrc.load_src('src_v2.pkl'); i=rng.integers(0,len(r),N)
    r=(r[i]-r[i].mean())*float(os.environ.get('WSCALE','1'))+R+(r[i].mean()-0.03); vr=vr[i]; vt=vt[i]; vz=vz[i]
    v87=np.sqrt(2*E*qe/(M87*amu)); v=np.sqrt(vr**2+vt**2+vz**2)
    # keep angles scaled by 1/sqrt(3) (same normalized emittance), speed for E
    ar=vr/vz*scale_ang; at=vt/vz*scale_ang
    sc=v87*np.sqrt(M87/M)
    vz=sc/np.sqrt(1+ar**2+at**2); vr=ar*vz; vt=at*vz
    qm=qe/(M*amu); az=2*np.pi*rng.random(N)
    bz0=np.interp(zs,s.zg,s.Bz[0]); vt=vt-qm*bz0*r/2
    ur=np.stack([np.cos(az),np.sin(az)],1); ut=np.stack([-np.sin(az),np.cos(az)],1)
    X=np.concatenate([r[:,None]*ur,np.full((N,1),zs)],1)
    V=np.concatenate([vr[:,None]*ur+vt[:,None]*ut,vz[:,None]],1)
    return X,V,np.full(N,qm)
def run_noSC(B0,N=3000):
    s=sol(B0); zp=np.arange(1.7,3.8,0.0025); out={}
    for k,(nm,M) in enumerate((('85',M85),('87',M87))):
        X,V,qm=beam(s,M,N,np.random.default_rng(10+k)); out[nm]=track(s,X,V,qm,zp,0.02/(qm*B0))
    k,r7,r5=analyze(out,zp); m=metrics(r7,r5); m.update(B0=B0,z=float(zp[k]),rms87=float(np.sqrt(np.mean(r7**2))),med85=float(np.median(r5))); return m
if __name__=='__main__':
    from multiprocessing import Pool
    sol(2.8)
    Bs=[float(x) for x in sys.argv[1:]]
    with Pool(4) as p: res=p.map(run_noSC,Bs)
    for r in res: print(r)
