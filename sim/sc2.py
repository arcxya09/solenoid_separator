# steady-state axisymmetric space charge: E_r(r,z)=lambda_enc(r,z)/(2 pi eps0 r), iterated
from ring import *
import ring, json, sys
from scan import metrics
import os
B0=float(os.environ.get("B0","2.74")); s=get_sol(B0)
zg=np.arange(-0.6,1.9,0.0025); rg=np.arange(0,0.09,0.0005)
def build_E(recs, ws):
    lam=np.zeros((len(zg),len(rg)))
    for rec,w in zip(recs,ws):
        r=np.hypot(rec[...,0],rec[...,1]); vz=rec[...,5]
        for i in range(len(zg)):
            ok=np.isfinite(r[i])
            idx=np.searchsorted(rg,r[i][ok]); idx=np.minimum(idx,len(rg)-1)
            np.add.at(lam[i],idx,w/vz[i][ok])
    lam=np.cumsum(lam,1)
    Er=np.zeros_like(lam); Er[:,1:]=lam[:,1:]/(2*np.pi*eps0*rg[1:])
    def Ef(X):
        r=np.hypot(X[:,0],X[:,1]); z=X[:,2]
        fz=np.clip((z-zg[0])/(zg[1]-zg[0]),0,len(zg)-1.001); fr=np.clip(r/(rg[1]-rg[0]),0,len(rg)-1.001)
        i=fz.astype(int); j=fr.astype(int); u=fz-i; t=fr-j
        e=Er[i,j]*(1-u)*(1-t)+Er[i+1,j]*u*(1-t)+Er[i,j+1]*(1-u)*t+Er[i+1,j+1]*u*t
        e=np.where((z<zg[0])|(z>zg[-1]),0,e); rs=np.where(r>0,r,1)
        return np.stack([e*X[:,0]/rs,e*X[:,1]/rs,0*e],1)
    return Ef
def run(I, N=2500, iters=4):
    ring.rng=np.random.default_rng(21)
    beams=[]; ab={'85':0.7217,'87':0.2783}
    for nm,M in (('85',M85),('87',M87)):
        X,V,qm=ring_beam(s,N,M,10e3); beams.append((nm,X,V,qm))
    Ef=None; zbest=None
    for it in range(iters if I>0 else 1):
        recs=[]; ws=[]
        for nm,X,V,qm in beams:
            recs.append(track(s,X,V,qm,zg,0.02/(qm*B0),Efield=Ef)); ws.append(I*ab[nm]/N)
        if I>0: Ef=build_E(recs,ws)
    out={'85':recs[0],'87':recs[1]}
    sel=zg>1.0
    r7a=np.hypot(recs[1][sel][...,0],recs[1][sel][...,1]); r5a=np.hypot(recs[0][sel][...,0],recs[0][sel][...,1])
    rms7=np.sqrt(np.nanmean(r7a**2,1)); k=np.nanargmin(rms7)
    m=metrics(r7a[k],r5a[k]); m.update(I=I, z=float(zg[sel][k]), rms87=float(rms7[k]))
    return m
if __name__=='__main__':
    from multiprocessing import Pool
    Is=[0,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4]
    with Pool(4) as p: res=p.map(run,Is)
    json.dump(res,open('sc.json','w'),indent=1)
    for r in res: print(r)
