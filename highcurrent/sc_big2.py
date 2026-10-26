import os,sys,json, numpy as np
from core import *
from scan import metrics
import big
B0=float(os.environ.get('B0','2.8')); s=big.sol(B0)
zg=np.arange(-0.6,3.5,0.0025); rg=np.arange(0,0.13,0.00025)
FN=float(os.environ.get('FN','0'))   # neutralization fraction (uniform)
def build_E(recs,ws):
    lam=np.zeros((len(zg),len(rg)))
    for rec,w in zip(recs,ws):
        r=np.hypot(rec[...,0],rec[...,1]); vz=rec[...,5]
        for i in range(len(zg)):
            ok=np.isfinite(r[i]); idx=np.minimum(np.searchsorted(rg,r[i][ok]),len(rg)-1)
            np.add.at(lam[i],idx,w/vz[i][ok])
    lam=np.cumsum(lam,1)*(1-FN)
    Er=np.zeros_like(lam); Er[:,1:]=lam[:,1:]/(2*np.pi*eps0*rg[1:])
    def Ef(X):
        r=np.hypot(X[:,0],X[:,1]); z=X[:,2]
        fz=np.clip((z-zg[0])/(zg[1]-zg[0]),0,len(zg)-1.001); fr=np.clip(r/(rg[1]-rg[0]),0,len(rg)-1.001)
        i=fz.astype(int); j=fr.astype(int); u=fz-i; t=fr-j
        e=Er[i,j]*(1-u)*(1-t)+Er[i+1,j]*u*(1-t)+Er[i,j+1]*(1-u)*t+Er[i+1,j+1]*u*t
        e=np.where((z<zg[0])|(z>zg[-1]),0,e); rs=np.where(r>0,r,1)
        return np.stack([e*X[:,0]/rs,e*X[:,1]/rs,0*e],1)
    return Ef
def run(I,N=2500,iters=6):
    ab={'85':0.7217,'87':0.2783}; beams=[]
    for k,(nm,M) in enumerate((('85',M85),('87',M87))):
        X,V,qm=big.beam(s,M,N,np.random.default_rng(30+k)); beams.append((nm,X,V,qm))
    Ef=None
    for it in range(iters if I>0 else 1):
        recs=[];ws=[]
        for nm,X,V,qm in beams:
            recs.append(track(s,X,V,qm,zg,0.02/(qm*B0),Efield=Ef)); ws.append(I*ab[nm]/N)
        if I>0:
            Efn=build_E(recs,ws)
            if Ef is None: Ef=Efn
            else:
                Eo=Ef; Ef=(lambda a,b: (lambda X: 0.5*a(X)+0.5*b(X)))(Eo,Efn)
    sel=zg>1.75
    r7a=np.hypot(recs[1][sel][...,0],recs[1][sel][...,1]); r5a=np.hypot(recs[0][sel][...,0],recs[0][sel][...,1])
    rms7=np.sqrt(np.nanmean(r7a**2,1)); k=np.nanargmin(rms7)
    m=metrics(r7a[k],r5a[k]); m.update(I=I,B0=B0,FN=FN,z=float(zg[sel][k]),rms87=float(rms7[k])); return m
if __name__=='__main__':
    print(json.dumps(run(float(sys.argv[1]))))
