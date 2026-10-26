# self-consistent axisymmetric space charge (sc2.build_E) with the real source distribution
from bl import *
import json, sys
from multiprocessing import Pool
zg=np.arange(-0.6,1.9,0.0025); rg=np.arange(0,0.09,0.0005)
def build_E(recs, ws):
    lam=np.zeros((len(zg),len(rg)))
    for rec,w in zip(recs,ws):
        r=np.hypot(rec[...,0],rec[...,1]); vz=rec[...,5]
        for i in range(len(zg)):
            ok=np.isfinite(r[i]); idx=np.minimum(np.searchsorted(rg,r[i][ok]),len(rg)-1)
            np.add.at(lam[i],idx,w/vz[i][ok])
    lam=np.cumsum(lam,1); Er=np.zeros_like(lam); Er[:,1:]=lam[:,1:]/(2*np.pi*eps0*rg[1:])
    def Ef(X):
        r=np.hypot(X[:,0],X[:,1]); z=X[:,2]
        fz=np.clip((z-zg[0])/(zg[1]-zg[0]),0,len(zg)-1.001); fr=np.clip(r/(rg[1]-rg[0]),0,len(rg)-1.001)
        i=fz.astype(int); j=fr.astype(int); u=fz-i; t=fr-j
        e=Er[i,j]*(1-u)*(1-t)+Er[i+1,j]*u*(1-t)+Er[i,j+1]*(1-u)*t+Er[i+1,j+1]*u*t
        e=np.where((z<zg[0])|(z>zg[-1]),0,e); rs=np.where(r>0,r,1)
        return np.stack([e*X[:,0]/rs,e*X[:,1]/rs,0*e],1)
    return Ef,Er
def run(a):
    I,B0,src,N=a
    s=get_sol(B0); beams={nm:beam(s,nm,N,src=src,seed=3) for nm in ('85','87')}
    Ef=None
    for it in range(5 if I>0 else 1):
        recs={nm:track(s,*beams[nm],zg,0.02/(beams[nm][2][0]*B0),Efield=Ef) for nm in beams}
        if I>0:
            Enew,Er=build_E([recs['85'],recs['87']],[I*AB['85']/N,I*AB['87']/N])
            Ef=Enew
    sel=zg>=1.0; out={nm:recs[nm][sel] for nm in recs}; zp=zg[sel]
    res=dict(I=I,B0=B0,src=src,best=best(out,zp))
    col=json.load(open('collector.json'))
    res['fixed']=fixed(out,col['z'],col['a'],zp)
    # T87 vs z for a fixed hole radius (focus shift)
    r7=np.hypot(out['87'][...,0],out['87'][...,1]); r5=np.hypot(out['85'][...,0],out['85'][...,1])
    res['zscan']=dict(z=zp[::4].round(4).tolist(),T=np.nanmean(r7<col['a'],1)[::4].round(4).tolist(),L=np.nanmean(r5<col['a'],1)[::4].round(5).tolist())
    if 'env' in sys.argv:
        r7a=np.hypot(recs['87'][...,0],recs['87'][...,1]); r5a=np.hypot(recs['85'][...,0],recs['85'][...,1])
        res['env']=dict(z=zg[::8].round(3).tolist(),r87=[np.nanpercentile(r7a[i],[5,50,95]).round(5).tolist() for i in range(0,len(zg),8)],
                        r85=[np.nanpercentile(r5a[i],[5,50,95]).round(5).tolist() for i in range(0,len(zg),8)])
    return res
if __name__=='__main__':
    cf=json.loads(open(sys.argv[1]).read())
    with Pool(4) as p: res=p.map(run,[tuple(c) for c in cf])
    json.dump(res,open(sys.argv[2],'w'))
    for r in res: print(r['I']*1e3,'mA',r['B0'],r['src'],'best',{k:round(v,4) for k,v in r['best'].items()},'fixed',{k:round(v,5) for k,v in r['fixed'].items()})
