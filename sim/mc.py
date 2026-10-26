from core import *
rng=np.random.default_rng(1)
def make_beam(N,M,E,zs=-0.6,rs=1e-3,th1=0.04,th2=0.08,sigE=2.0,sig_th_therm=0.0):
    m,v0,qm=ion(M,E)
    r=rs*np.sqrt(rng.random(N)); ph=2*np.pi*rng.random(N)
    th=np.sqrt(th1**2+(th2**2-th1**2)*rng.random(N))   # uniform in solid angle
    az=2*np.pi*rng.random(N)
    Ek=E+sigE*rng.standard_normal(N); v=np.sqrt(2*Ek*qe/m)
    X=np.stack([r*np.cos(ph),r*np.sin(ph),np.full(N,zs)],1)
    d=np.stack([np.sin(th)*np.cos(az),np.sin(th)*np.sin(az),np.cos(th)],1)
    if sig_th_therm>0:
        d[:,:2]+=sig_th_therm*rng.standard_normal((N,2)); d/=np.linalg.norm(d,axis=1)[:,None]
    return X,d*v[:,None],np.full(N,qm),th
def run(sol,N=4000,E=10e3,zplanes=np.arange(1.0,1.8,0.005),**kw):
    out={}
    for nm,M in (('85',M85),('87',M87)):
        X,V,qm,th=make_beam(N,M,E,**kw)
        dt=0.02/(qm*sol.B0)
        rec=track(sol,X,V,qm,zplanes,dt)
        out[nm]=(rec,th)
    return out
if __name__=='__main__':
    import time; t=time.time()
    s=get_sol(2.9); zp=np.arange(1.0,1.8,0.005)
    out=run(s,N=3000,zplanes=zp)
    print('time',time.time()-t)
    r7=np.hypot(out['87'][0][...,0],out['87'][0][...,1]); r5=np.hypot(out['85'][0][...,0],out['85'][0][...,1])
    rms7=np.sqrt(np.nanmean(r7**2,1)); 
    k=np.nanargmin(rms7); print('best z',zp[k],'rms87 mm',rms7[k]*1e3,'r85 min/median',np.nanmin(r5[k])*1e3,np.nanmedian(r5[k])*1e3)
    for q in [50,90,99,99.9]: print(q, np.nanpercentile(r7[k],q)*1e3)
    np.savez('mc_base.npz',zp=zp,r7=r7,r5=r5)
