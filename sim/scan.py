from ring import *
import json,sys
from multiprocessing import Pool
def metrics(r7,r5):
    ra,n7,n5,pur=purity_curve(r7,r5)
    res={}
    for P in (0.99,0.999):
        ok=pur>=P; res['T_%g'%P]=float(n7[ok].max()) if ok.any() else 0.0
        res['a_%g'%P]=float(ra[ok][np.argmax(n7[ok])]) if ok.any() else 0.0
    return res
def job(cfg):
    global rng
    import ring; ring.rng=np.random.default_rng(cfg.get('seed',3))
    B0=cfg.pop('B0'); E=cfg.pop('E',10e3); N=cfg.pop('N',3000); zp=np.arange(1.0,2.2,0.0025)
    s=get_sol(B0); out=run_ring(s,N=N,E=E,zplanes=zp,**{k:v for k,v in cfg.items() if k not in ('seed','tag')})
    k,r7,r5=analyze(out,zp)
    m=metrics(r7,r5); m.update(z=float(zp[k]),B0=B0,E=E,med85=float(np.median(r5)),r7_90=float(np.percentile(r7,90)))
    m.update({k:v for k,v in cfg.items()})
    return m
if __name__=='__main__':
    cfgs=json.loads(open(sys.argv[1]).read())
    with Pool(4) as p: res=p.map(job,cfgs)
    json.dump(res,open(sys.argv[2],'w'),indent=1)
    for r in res: print(r)
