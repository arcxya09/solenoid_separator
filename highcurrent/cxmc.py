# charge-exchange fast neutrals: where do they land? (big design, B=2.8 T)
import numpy as np, json
from core import *
import big
s=big.sol(2.8); zf=2.335; a_hole=5.35e-3
zp=np.arange(-0.595,zf,0.005); res={}
for k,(nm,M) in enumerate((('85',M85),('87',M87))):
    X,V,qm=big.beam(s,M,2000,np.random.default_rng(40+k))
    rec=track(s,X,V,qm,zp,0.02/(qm*2.8))
    x,y,z,vx,vy,vz=[rec[...,i] for i in range(6)]
    t=(zf-z)/vz; xf=x+vx*t; yf=y+vy*t; rf=np.hypot(xf,yf)
    dl=np.sqrt(vx**2+vy**2+vz**2)/vz*0.005           # path length per plane step
    hit=(rf<a_hole)
    w=dl/np.nansum(dl,0,keepdims=True)                 # fraction of path per element
    # fraction of CX events (uniform per path length) whose neutral lands in hole, as function of z
    prof=np.nansum(w*hit,1)/rec.shape[1]
    res[nm]=dict(p_hit_total=float(np.nansum(prof)), prof=prof.tolist(), Lpath=float(np.nanmean(np.nansum(dl,0))))
    print(nm,'p_hit',res[nm]['p_hit_total'],'path',res[nm]['Lpath'])
res['z']=zp.tolist()
# purity vs pressure
kB=1.380649e-23
out=[]
for sig in (1e-16,3e-16,1e-15):
  for P_torr in (2e-6,5e-6,1e-5,2e-5,5e-5):
    n=P_torr*133.322/(kB*300)*1e-6     # cm^-3
    f5=1-np.exp(-n*sig*res['85']['Lpath']*100); f7=1-np.exp(-n*sig*res['87']['Lpath']*100)
    imp=0.7217*f5*res['85']['p_hit_total']; good=0.2783*(1-f7)*1.0 + 0.2783*f7*res['87']['p_hit_total']
    out.append(dict(sig=sig,P=P_torr,f_cx=f5,purity=good/(good+imp),loss87=f7*(1-res['87']['p_hit_total'])))
res['pur']=out
json.dump(res,open('cx.json','w'))
for o in out: print(o)
