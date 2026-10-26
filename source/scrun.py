from es import *
import json,sys,pickle
from multiprocessing import Pool
BASE=dict(flare=67.5,h=1e-3,zl=30e-3,zhi=90e-3)
def run(I,Vl=-15000,N=2000,iters=6,seed=1,save=None,extra=None):
    c=dict(BASE,Vl=Vl); c.update(extra or {})
    g=Geom(**c); g.solve()
    X,V,m=emit(g,N,rng=np.random.default_rng(seed))
    rho=None
    for it in range(iters if I>0 else 1):
        out,alive,rnew=track(g,X,V,m,zout=0.085,I=I,deposit=I>0)
        if I>0:
            rho=rnew if rho is None else 0.5*rho+0.5*rnew
            g.solve(rho)
    out,alive,_=track(g,X,V,m,zout=0.085)
    r,_=beam_stats(out,g); r={k:float(v) for k,v in r.items()}
    r.update(I=I,Vl=Vl,trans=float(np.isfinite(out[:,0]).mean()))
    ok=np.isfinite(out[:,0]); x,y,z,vx,vy,vz=out[ok].T; rr=np.hypot(x,y)
    r['tp_rms']=float(np.std((x*vy-y*vx)/rr/vz)); r['E_rms']=float(np.std(0.5*m*(vx**2+vy**2+vz**2)/qe))
    if save:
        pickle.dump(dict(out=out,phi=g.phi[::4,::4],r=g.r[::4],z=g.z[::4],solid=g.solid[::4,::4],rho=rho,cfg=c,I=I,m=m),open(save,'wb'))
        # trajectories for plotting
    return r
def job(a): return run(*a)
if __name__=='__main__':
    args=json.loads(sys.argv[1])
    with Pool(4) as p: res=p.map(job,[tuple(a) for a in args])
    json.dump(res,open(sys.argv[2],'w'))
    for r in res: print(r['I']*1e3,'mA Vl',r['Vl'],'trans',round(r['trans'],3),'slope',round(r['slope'],2),'rp_mean',round(r['rp_mean']*1e3,2),'rp_rms',round(r['rp_rms']*1e3,2),'tp',round(r['tp_rms']*1e3,2),'wid',round(r['rwid_rms']*1e3,2),'eps',round(r['eps_rms']*1e6,2))
