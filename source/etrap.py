# electrons born (residual-gas ionization) inside the +V einzel of the v3 source: are they trapped by the residual solenoid field?
from es import *
import json,sys,os
from multiprocessing import Pool
me=9.1093837e-31; QM=-qe/me
V3=dict(flare=67.5,h=1e-3,zl=30e-3,zhi=90e-3,rnd=0.0014,t=0.003,d=0.010,Vmid=-2500,d2=0.006,t2=0.003,sl2=0.004,sl=0.006,gl=0.010)
def run(a):
    B0,dlnB,Vl,N,T,extra=a
    g=Geom(**dict(V3,Vl=Vl,**extra)); g.solve(); rng=np.random.default_rng(1)
    zl=g.zl; zc=zl+g.Ll/2
    X=[];
    while sum(len(x) for x in X)<N:
        r=g.R+0.0006*rng.standard_normal(4*N); z=rng.uniform(zl-g.gl,zl+g.Ll+g.gl,4*N)
        ph,_=g.interp(g.phi,r,z); sm,_=g.interp(g.solidf,r,z); ok=(ph>20)&(sm<0.01)
        X.append(np.stack([r[ok],0*r[ok],z[ok]],1))
    X=np.concatenate(X)[:N]; ph0,_=g.interp(g.phi,X[:,0],X[:,2])
    v=np.sqrt(2*3*qe/me); u=rng.standard_normal((N,3)); V=v*u/np.linalg.norm(u,axis=1)[:,None]
    alive=np.ones(N,bool); fate=np.zeros(N,int); tloss=np.full(N,np.nan); dt=1e-12
    ts=[];surv=[]
    for n in range(int(T/dt)):
        a=np.where(alive)[0]
        if len(a)==0: break
        x,y,z=X[a].T; r=np.hypot(x,y)
        Er,_=g.interp(g.Er,r,z); Ez,_=g.interp(g.Ez,r,z)
        E=np.stack([Er*x/r,Er*y/r,Ez],1)
        Bz=B0*(1+dlnB*(z-zc)); Br=-0.5*r*B0*dlnB
        B=np.stack([Br*x/r,Br*y/r,Bz],1)
        vm=V[a]+0.5*QM*dt*E; tt=0.5*QM*dt*B; s2=2*tt/(1+np.sum(tt*tt,1))[:,None]
        vp=vm+np.cross(vm,tt); vm=vm+np.cross(vp,s2); Vn=vm+0.5*QM*dt*E
        Xn=X[a]+Vn*dt; X[a]=Xn; V[a]=Vn
        rr=np.hypot(Xn[:,0],Xn[:,1]); sm,_=g.interp(g.solidf,rr,Xn[:,2])
        hit=(sm>0.5)|(rr<g.r[0])|(rr>g.r[-1])|(Xn[:,2]<g.z[0])|(Xn[:,2]>g.z[-1])
        if hit.any():
            ids=a[hit]; alive[ids]=False; tloss[ids]=n*dt
            phh,_=g.interp(g.phi,np.clip(rr[hit],g.r[0],g.r[-1]),np.clip(Xn[hit,2],g.z[0],g.z[-1])); fate[ids]=np.where(phh>Vl/2,1,2)
        if n%5000==0: ts.append(n*dt); surv.append(float(alive.mean()))
    return dict(B0=B0,dlnB=dlnB,Vl=Vl,extra=extra,T=T,surv_end=float(alive.mean()),to_einzel=float((fate==1).mean()),to_other=float((fate==2).mean()),
                t_med=float(np.nanmedian(tloss)) if np.isfinite(tloss).any() else None,ts=ts,surv=surv,phi_birth_med=float(np.median(ph0)))
if __name__=='__main__':
    args=json.loads(sys.argv[1])
    with Pool(int(os.environ.get('NP','4'))) as p:
        for r in p.imap_unordered(run,[tuple(a) for a in args]):
            print({k:v for k,v in r.items() if k not in('ts','surv')},flush=True)
            try: old=json.load(open(sys.argv[2]))
            except Exception: old=[]
            old.append(r); json.dump(old,open(sys.argv[2],'w'))
