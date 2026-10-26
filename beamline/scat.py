# charge exchange (fast neutrals, straight line) and small-angle elastic scattering on residual gas:
# probability that an ion suffering one event at position z lands inside the Rb-87 hole
from bl import *
import json, sys
from multiprocessing import Pool
col=json.load(open('collector.json')); ZC=col['z']; A=col['a']; B0=col.get('B0',2.78)
ZE=np.arange(-0.55,ZC-0.02,0.05)            # event planes
ANG=[0.0,2e-3,5e-3,10e-3,20e-3,50e-3]       # 0 = charge exchange (neutral)
def job(a):
    nm,j,ang=a; s=get_sol(B0); N=4000
    X,V,qm=beam(s,nm,N,seed=11)
    zp=np.array([ZE[j],ZC])
    if ZE[j]>-0.6+1e-6:
        rec=track(s,X,V,qm,zp[:1],0.02/(qm[0]*B0)); P=rec[0]; X=P[:,:3].copy(); V=P[:,3:].copy()
    if ang==0.0:   # neutral: straight line
        t=(ZC-X[:,2])/V[:,2]; xf=X[:,0]+V[:,0]*t; yf=X[:,1]+V[:,1]*t
    else:
        rng=np.random.default_rng(j); v=np.linalg.norm(V,axis=1); u=V/v[:,None]
        e1=np.cross(u,[0,0,1.]); e1/=np.linalg.norm(e1,axis=1)[:,None]; e2=np.cross(u,e1)
        ph=2*np.pi*rng.random(N)
        V=v[:,None]*(np.cos(ang)*u+np.sin(ang)*(np.cos(ph)[:,None]*e1+np.sin(ph)[:,None]*e2))
        rec=track(s,X,V,qm,np.array([ZC]),0.02/(qm[0]*B0)); xf,yf=rec[0][:,0],rec[0][:,1]
    r=np.hypot(xf,yf); return dict(nm=nm,z=float(ZE[j]),ang=ang,p_hole=float(np.nanmean(r<A)),p_lost=float(np.mean(~np.isfinite(r))))
if __name__=='__main__':
    cf=[(nm,j,a) for nm in ('85','87') for j in range(len(ZE)) for a in ANG]
    with Pool(4) as p: res=p.map(job,cf)
    json.dump(dict(ZC=ZC,A=A,res=res),open('scat.json','w'))
    for a in ANG:
        for nm in ('85','87'):
            v=[r['p_hole'] for r in res if r['nm']==nm and r['ang']==a]; print(nm,a,'mean p_hole %.4f'%np.mean(v),[round(x,3) for x in v])
