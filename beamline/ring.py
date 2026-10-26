from core import *
rng=np.random.default_rng(7)
def ring_beam(sol,N,M,E,R=0.03,w=1e-3,sigp=4e-3,sigE=1.0,zs=-0.6,Bs=0.0,naz=None,az0=0.0):
    m,v0,qm=ion(M,E)
    r=R+w*(rng.random(N)-0.5)
    az=2*np.pi*rng.random(N) if naz is None else np.full(N,az0)
    Ek=E+sigE*rng.standard_normal(N); v=np.sqrt(2*Ek*qe/m)
    xp=sigp*rng.standard_normal(N); yp=sigp*rng.standard_normal(N)  # radial / azimuthal angle (local frame)
    # zero canonical angular momentum (source in field-free region) + residual Bs
    bz0=np.interp(zs,sol.zg,sol.Bz[0])
    vth_add=qm*(Bs*r/2 - bz0*r/2)
    ur=np.stack([np.cos(az),np.sin(az)],1); ut=np.stack([-np.sin(az),np.cos(az)],1)
    vt=v[:,None]*(xp[:,None]*ur+yp[:,None]*ut)+vth_add[:,None]*ut
    vz=np.sqrt(v**2-np.sum(vt**2,1))
    X=np.stack([r*np.cos(az),r*np.sin(az),np.full(N,zs)],1)
    V=np.concatenate([vt,vz[:,None]],1)
    return X,V,np.full(N,qm)
def run_ring(sol,N=3000,E=10e3,zplanes=np.arange(1.0,1.6,0.0025),**kw):
    out={}
    for nm,M in (('85',M85),('87',M87)):
        X,V,qm=ring_beam(sol,N,M,E,**kw)
        out[nm]=track(sol,X,V,qm,zplanes,0.02/(qm*sol.B0))
    return out
def analyze(out,zp,which='87'):
    r7=np.hypot(out['87'][...,0],out['87'][...,1]); r5=np.hypot(out['85'][...,0],out['85'][...,1])
    tgt=r7 if which=='87' else r5
    rms=np.sqrt(np.nanmean(tgt**2,1)); k=np.nanargmin(rms)
    return k,r7[k],r5[k]
def purity_curve(r7,r5,ab=0.2783):
    ra=np.linspace(0,0.02,401)
    n7=np.array([(r7<a).mean() for a in ra]); n5=np.array([(r5<a).mean() for a in ra])
    pur=ab*n7/np.maximum(ab*n7+(1-ab)*n5,1e-30)
    return ra,n7,n5,pur
if __name__=='__main__':
    import time,sys; t=time.time()
    B0=float(sys.argv[1]) if len(sys.argv)>1 else 2.85
    s=get_sol(B0); zp=np.arange(1.0,1.6,0.0025)
    out=run_ring(s,N=3000,zplanes=zp)
    k,r7,r5=analyze(out,zp); print('t',time.time()-t,'zbest',zp[k])
    print('r87 pct 50/90/99/99.9 mm',[round(np.percentile(r7,q)*1e3,2) for q in (50,90,99,99.9)])
    print('r85 pct 0.1/1/10/50 mm',[round(np.percentile(r5,q)*1e3,2) for q in (0.1,1,10,50)])
    ra,n7,n5,pur=purity_curve(r7,r5)
    for a in [2,3,4,5,6]:
        i=np.argmin(abs(ra-a*1e-3)); print(a,'mm: T87',n7[i],'leak85',n5[i],'purity',pur[i])
