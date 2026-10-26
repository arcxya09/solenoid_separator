# machine B, 6 mA unneutralised: for each B find best collector plane (free z) -> choose design z, B, hole
from hc30 import *
import sys
def job2(a):
    B,I,src,N=a; wp.SRC=src; nm='h6_%g_%g'%(B,I); mk30(nm,B)
    zp=np.arange(ZC-0.40,ZC+0.15,0.005)
    recs,zp,_=run_beam(nm,src,N=N,I=I*1e-3,zplanes=zp,seed=1,beam_kw=dict(dR=0.0),iters=5)
    out=[]
    for k in range(len(zp)):
        b=at_plane(recs,zp,k,0.999,recenter=True,amax=0.025); out.append((float(zp[k]),b['T'],b['a'],b['rms87'],b['r85min'],b['r87max']))
    np.save('recs_%s.npy'%nm,{kk:v for kk,v in recs.items()},allow_pickle=True)
    return dict(B=B,I=I,src=src,planes=out)
if __name__=='__main__':
    src=sys.argv[1]; I=float(sys.argv[2]); Bs=[float(x) for x in sys.argv[3].split(',')]
    with Pool(4) as p: res=p.map(job2,[(B,I,src,2000) for B in Bs],chunksize=1)
    json.dump(res,open('hc6_%g.json'%I,'w'))
    for r in res:
        pl=[x for x in r['planes'] if x[1]>=0.99]
        if pl:
            zs=[x[0] for x in pl]; best=max(r['planes'],key=lambda x:(x[1],-x[2]))
            print('I=%g B=%.4f  T>=99%% for z in [%.3f,%.3f] (%d planes); best z=%.3f a=%.1f mm rms=%.1f r85min=%.1f r87max=%.1f'%(r['I'],r['B'],min(zs),max(zs),len(zs),best[0],best[2]*1e3,best[3]*1e3,best[4]*1e3,best[5]*1e3),flush=True)
        else:
            best=max(r['planes'],key=lambda x:x[1]); print('I=%g B=%.4f FAIL best T=%.3f at z=%.3f'%(r['I'],r['B'],best[1],best[0]),flush=True)
