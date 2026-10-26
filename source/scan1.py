from es import *
import json,sys,itertools
from multiprocessing import Pool
def job(c):
    g=Geom(**c); g.solve(); X,V,m=emit(g,600,rng=np.random.default_rng(1))
    out,alive,_=track(g,X,V,m,zout=0.065); r,_=beam_stats(out,g)
    r={k:float(v) for k,v in r.items()}; r.update(c); r['trans']=float(np.isfinite(out[:,0]).mean()); return r
cfgs=[dict(rc=rc,flare=fl,h=h) for rc in (None,3e-3,1.5e-3,0.8e-3) for fl in (0,67.5) for h in (0.2e-3,1.0e-3)]
with Pool(4) as p: res=p.map(job,cfgs)
json.dump(res,open('scan1.json','w'))
for r in res: print(r['rc'],r['flare'],r['h'],'trans',round(r['trans'],2),'slope',round(r['slope'],1),'wid',round(r['rwid_rms']*1e3,2),'unc',round(r['rp_uncorr']*1e3,2),'eps',round(r['eps_rms']*1e6,2))
