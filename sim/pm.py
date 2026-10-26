from ring import *
import ring, json, sys
from scan import metrics
from multiprocessing import Pool
KW=dict(L=2.0,a_in=0.09,a_out=0.12,nz=400,nr=3,rmax=0.085,zmin=-1.0,zmax=4.5)
def job(B0):
    ring.rng=np.random.default_rng(3)
    s=get_sol(B0,**KW); zp=np.arange(2.0,3.6,0.0025)
    out=run_ring(s,N=2500,E=3.3e3,zplanes=zp,R=0.045,sigp=7e-3,w=1e-3,sigE=1.0,zs=-0.6)
    k,r7,r5=analyze(out,zp); m=metrics(r7,r5); m.update(B0=B0,z=float(zp[k]),med85=float(np.median(r5)),r7_90=float(np.percentile(r7,90)))
    return m
if __name__=='__main__':
    get_sol(1.2,**KW)
    Bs=[float(x) for x in sys.argv[1:]]
    with Pool(4) as p: res=p.map(job,Bs)
    for r in res: print(r)
    json.dump(res,open('pm_%s.json'%sys.argv[1],'w'))
