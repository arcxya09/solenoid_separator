# feed the simulated ring-source phase space into the solenoid separator
from ring import *
import ring, pickle, json, sys
from scan import metrics
from multiprocessing import Pool
SRC='/tmp/claude-0/-home-claude/f9c16c93-b28c-5a0f-9f67-9ad570622977/scratchpad/src/'
def load_src(fn):
    d=pickle.load(open(SRC+fn,'rb')); o=d['out']; ok=np.isfinite(o[:,0]); o=o[ok]
    x,y,z,vx,vy,vz=o.T; r=np.hypot(x,y)
    vr=(x*vx+y*vy)/r; vt=(x*vy-y*vx)/r
    return r,vr,vt,vz
def src_beam(sol,M,srcfile,N,zs=-0.6,dx=0.0,dy=0.0,tilt=0.0,dR=0.0,rng=None):
    rng=rng or np.random.default_rng(0)
    r,vr,vt,vz=load_src(srcfile); i=rng.integers(0,len(r),N)
    r,vr,vt,vz=r[i]+dR,vr[i],vt[i],vz[i]
    # isotope: source simulated for Rb-87; same voltage => velocities scale sqrt(m87/m)
    sc=np.sqrt(M87/M); vr,vt,vz=vr*sc,vt*sc,vz*sc
    qm=qe/(M*amu)
    az=2*np.pi*rng.random(N)
    bz0=np.interp(zs,sol.zg,sol.Bz[0]); vt=vt-qm*bz0*r/2     # source region shielded: zero canonical momentum
    ur=np.stack([np.cos(az),np.sin(az)],1); ut=np.stack([-np.sin(az),np.cos(az)],1)
    X=np.concatenate([r[:,None]*ur+np.array([dx,dy]),np.full((N,1),zs)],1)
    vxy=vr[:,None]*ur+vt[:,None]*ut
    V=np.concatenate([vxy,vz[:,None]],1)
    if tilt:  # tilt the whole source about y axis
        c,s=np.cos(tilt),np.sin(tilt); V=np.stack([c*V[:,0]+s*V[:,2],V[:,1],-s*V[:,0]+c*V[:,2]],1)
    return X,V,np.full(N,qm)
def job(cfg):
    cfg=dict(cfg); B0=cfg.pop('B0'); fn=cfg.pop('src'); N=cfg.pop('N',3000)
    s=get_sol(B0); zp=np.arange(1.0,2.2,0.0025); out={}
    for k,(nm,M) in enumerate((('85',M85),('87',M87))):
        X,V,qm=src_beam(s,M,fn,N,rng=np.random.default_rng(10+k),**cfg)
        out[nm]=track(s,X,V,qm,zp,0.02/(qm*B0))
    r7a=np.hypot(out['87'][...,0],out['87'][...,1]); r5a=np.hypot(out['85'][...,0],out['85'][...,1])
    # allow the collector hole to be re-centred on the Rb-87 spot centroid (handles source offset)
    rms7=np.sqrt(np.nanmean(r7a**2,1)); k=np.nanargmin(rms7)
    cx=np.nanmean(out['87'][k][:,0]); cy=np.nanmean(out['87'][k][:,1])
    r7=np.hypot(out['87'][k][:,0]-cx,out['87'][k][:,1]-cy); r5=np.hypot(out['85'][k][:,0]-cx,out['85'][k][:,1]-cy)
    m=metrics(r7,r5); m.update(B0=B0,src=fn,z=float(zp[k]),cx=float(cx),cy=float(cy),rms87=float(np.std(r7)),med85=float(np.median(r5)),**cfg)
    if 'save' in sys.argv: m['xy']={nm:[(out[nm][k][:2500,0]*1e3).round(3).tolist(),(out[nm][k][:2500,1]*1e3).round(3).tolist()] for nm in out}
    return m
if __name__=='__main__':
    cfgs=json.loads(sys.argv[1])
    with Pool(4) as p: res=p.map(job,cfgs)
    json.dump(res,open(sys.argv[2],'w'))
    for r in res: print({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items() if k!='xy'})
