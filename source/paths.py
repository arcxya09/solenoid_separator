# trajectories for plotting, using saved self-consistent 0.1 mA charge density
from es import *
import pickle,json
d=pickle.load(open('src_0p1mA.pkl','rb'))
g=Geom(**d['cfg']); 
# rebuild full-res rho not saved at full res -> recompute quickly with 2 iterations
X,V,m=emit(g,2000,rng=np.random.default_rng(4)); g.solve()
rho=None
for it in range(4):
    out,al,rn=track(g,X,V,m,0.085,I=1e-4,deposit=True); rho=rn if rho is None else 0.5*rho+0.5*rn; g.solve(rho)
# record paths of 40 particles
Xp,Vp,_=emit(g,40,rng=np.random.default_rng(5)); qm=qe/m
paths=[[] for _ in range(40)]
Xc=Xp.copy(); Vc=Vp.copy()
for n in range(200000):
    r=np.hypot(Xc[:,0],Xc[:,1]); Er,_=g.interp(g.Er,r,Xc[:,2]); Ez,_=g.interp(g.Ez,r,Xc[:,2])
    vm=np.linalg.norm(Vc,axis=1); dt=np.clip(0.02*g.dr/np.maximum(vm,1),2e-11,5e-9)
    Vc=Vc+qm*np.stack([Er*Xc[:,0]/r,Er*Xc[:,1]/r,Ez],1)*dt[:,None]; Xc=Xc+Vc*dt[:,None]
    if n%40==0:
        for i in range(40): paths[i].append([float(Xc[i,2]*1e3),float((np.hypot(Xc[i,0],Xc[i,1])-g.R)*1e3)])
    if (Xc[:,2]>0.088).all(): break
# thin paths
P=[]
for p in paths:
    p=np.array(p); keep=[0]
    for k in range(1,len(p)):
        if np.hypot(*(p[k]-p[keep[-1]]))>0.25: keep.append(k)
    P.append(p[keep].round(3).tolist())
sub=(slice(None,None,2),slice(None,None,4))
dd=dict(z=(g.z[::4]*1e3).round(2).tolist(), x=((g.r[::2]-g.R)*1e3).round(2).tolist(),
        phi=(g.phi[sub]/1e3).round(3).T.tolist(), solid=g.solid[sub].astype(int).T.tolist(), paths=P)
json.dump(dd,open('paths.json','w')); print(len(dd['z']),len(dd['x']))
