from es import *
import json,sys
from multiprocessing import Pool
B=dict(flare=67.5,h=1e-3,zl=30e-3,zhi=90e-3)
def job(c):
    cc=dict(B); cc.update(c)
    g=Geom(**cc); g.solve(); X,V,m=emit(g,600,rng=np.random.default_rng(1))
    out,alive,_=track(g,X,V,m,zout=0.085); r,_=beam_stats(out,g)
    r={k:float(v) for k,v in r.items()}; r.update(c); r['trans']=float(np.isfinite(out[:,0]).mean()); return r
if __name__=='__main__':
    cfgs=[dict(zip(('Vl','loff'),v)) for v in json.loads(sys.argv[1])]
    with Pool(4) as p: res=p.map(job,cfgs)
    for r in res: print(r['Vl'],r.get('loff'),'trans',round(r['trans'],2),'slope',round(r['slope'],2),'rp_mean',round(r['rp_mean']*1e3,2),'rp_rms',round(r['rp_rms']*1e3,2),'wid',round(r['rwid_rms']*1e3,2),'unc',round(r['rp_uncorr']*1e3,2),'eps',round(r['eps_rms']*1e6,2))
