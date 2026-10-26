from core import *
s=get_sol(3.0)
E=10e3; zs=-0.6; th=0.06
zpl=np.arange(-0.55,3.0,0.01)
for B0 in [2.6,2.8,3.0,3.2,3.4]:
    s.set_B(B0)
    out=[]
    for M in (M85,M87):
        m,v,qm=ion(M,E)
        X0=np.array([[0,0,zs]]); V0=v*np.array([[np.sin(th),0,np.cos(th)]])
        rec=track(s,X0,V0,np.array([qm]),zpl,dt=0.01/(qm*B0))
        r=np.hypot(rec[:,0,0],rec[:,0,1]); out.append(r)
    # find axis crossings after z>1.0 (minima)
    for nm,r in zip(('85','87'),out):
        idx=[i for i in range(1,len(r)-1) if r[i]<r[i-1] and r[i]<=r[i+1] and zpl[i]>-0.3]
        print(B0,nm,[(round(zpl[i],2),round(r[i]*1e3,2)) for i in idx][-4:], 'rmax',round(np.nanmax(r)*1e3,1))
