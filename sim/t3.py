from mc import *
s=get_sol(2.9); zp=np.arange(1.0,1.8,0.005)
for th in [0.02,0.04,0.06,0.08]:
  for nm,M in (('87',M87),('85',M85)):
    m,v,qm=ion(M,10e3)
    X=np.array([[0,0,-0.6]]); V=v*np.array([[np.sin(th),0,np.cos(th)]])
    rec=track(s,X,V,np.array([qm]),zp,0.02/(qm*2.9))
    r=np.hypot(rec[:,0,0],rec[:,0,1]); k=np.nanargmin(r)
    print(th,nm,'min r at',round(zp[k],3),round(r[k]*1e3,3),' r@1.437',round(r[np.argmin(abs(zp-1.437))]*1e3,2))
