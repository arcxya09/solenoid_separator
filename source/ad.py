# accel-decel source scan: emitter +10 kV, suppressor Vmid<0, ground, optional einzel (either polarity)
from scrun import *
from hv import surf
import json,sys
from multiprocessing import Pool
def job(c):
    c=dict(c); I=c.pop('I',1e-4)
    full=dict(BASE,**c); g=Geom(**full); g.solve(); sf=surf(g)
    pl=g.interp(g.phi,np.full(len(g.z),g.R),g.z)[0]; zz=g.z
    after=zz>g.dg+g.t
    r=run(I,full.get('Vl'),N=1500,seed=3,extra={k:v for k,v in c.items() if k!='Vl'})
    sup=(zz>=g.d)&(zz<=g.dg+g.t)
    return dict(cfg=c,surf=sf,phisup=float(pl[sup].min()),phimin=float(pl[after].min()),phimax=float(pl[after].max()),
                **{k:r[k] for k in ('slope','rp_rms','rp_mean','tp_rms','rwid_rms','eps_rms','trans')})
def show(r):
    print(r['cfg'],'Emax',{k:round(v,2) for k,v in r['surf'].items()},'phi[%d,%d] sup %d'%(r['phimin'],r['phimax'],r.get('phisup',0)),'slope',round(r['slope'],2),'rp',round(r['rp_rms']*1e3,2),'mean',round(r['rp_mean']*1e3,2),'wid',round(r['rwid_rms']*1e3,2),'eps',round(r['eps_rms']*1e6,2),'tr',round(r['trans'],3),flush=True)
if __name__=='__main__':
    args=json.loads(sys.argv[1]); out=sys.argv[2]
    with Pool(int(sys.argv[3]) if len(sys.argv)>3 else 2) as p:
        for r in p.imap_unordered(job,args):
            show(r)
            try: old=json.load(open(out))
            except Exception: old=[]
            old.append(r); json.dump(old,open(out,'w'))
