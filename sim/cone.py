# original concept: point source on axis + oblique cone; focal position vs entry angle
from core import *
import json
s=get_sol(2.9); zp=np.arange(0.9,1.8,0.0025)
ths=np.linspace(0.01,0.08,15); res={}
for nm,M in (('85',M85),('87',M87)):
    m,v,qm=ion(M,10e3); N=len(ths)
    X=np.zeros((N,3)); X[:,2]=-0.6
    V=v*np.stack([np.sin(ths),0*ths,np.cos(ths)],1)
    rec=track(s,X,V,np.full(N,qm),zp,0.02/(qm*2.9))
    r=np.hypot(rec[...,0],rec[...,1])
    zf=[]
    for j in range(N):
        # last axis crossing (min of r) after z>0.95
        rr=r[:,j]; i=np.nanargmin(rr); zf.append(float(zp[i]))
    res[nm]=zf
res['theta']=ths.tolist()
json.dump(res,open('cone.json','w')); print(res)
