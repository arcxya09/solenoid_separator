# peak surface field vs edge rounding, beam quality check, potential barrier on beam line
from scrun import *
import json,sys
from multiprocessing import Pool
def surf(g):
    E=np.hypot(g.Er,g.Ez); s=g.solid
    nb=np.zeros_like(s)
    for a,b in ((1,0),(-1,0),(0,1),(0,-1)): nb|=np.roll(s,(a,b),(0,1))
    edge=(~g.fix)&nb
    edge[:,:3]=False; edge[:,-3:]=False; edge[:3,:]=False; edge[-3:,:]=False
    res={}
    for m,v,nm in g.parts:
        mm=np.zeros_like(s)
        for a,b in ((1,0),(-1,0),(0,1),(0,-1)): mm|=np.roll(m,(a,b),(0,1))
        sel=edge&mm
        res[nm]=float(E[sel].max()/1e6) if sel.any() else 0.0
    return res
def job(a):
    rnd,Vl,ex=a
    c=dict(BASE,Vl=Vl,rnd=rnd,**ex); g=Geom(**c); g.solve()
    sf=surf(g)
    # beam-line potential minimum (electron barrier)
    phi_line=g.interp(g.phi,np.full(len(g.z),g.R),g.z)[0]
    r=run(1e-4,Vl,N=1500,seed=3,extra=dict(rnd=rnd,**ex))
    return dict(rnd=rnd,Vl=Vl,ex=ex,surf=sf,phimin=float(phi_line[g.z>0.025].min()),rp_rms=r['rp_rms'],rp_mean=r['rp_mean'],tp=r['tp_rms'],wid=r['rwid_rms'],eps=r['eps_rms'],trans=r['trans'])
if __name__=='__main__':
    args=json.loads(sys.argv[1])
    with Pool(4) as p: res=p.map(job,[tuple(a) for a in args])
    json.dump(res,open(sys.argv[2],'w'))
    for r in res: print(r['rnd']*1e3,r['Vl'],r['ex'],{k:round(v,2) for k,v in r['surf'].items()},'phimin',round(r['phimin']),'rp',round(r['rp_rms']*1e3,2),'mean',round(r['rp_mean']*1e3,2),'wid',round(r['wid']*1e3,2),'eps',round(r['eps']*1e6,2),'tr',r['trans'])
