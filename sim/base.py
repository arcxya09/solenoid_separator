from ring import *
import ring, json
B0=2.74; s=get_sol(B0); out={}
z,b=s.onaxis(); sel=(z>-0.8)&(z<2.2)
out['field']=dict(z=z[sel][::5].round(4).tolist(), b=b[sel][::5].round(5).tolist())
# trajectories (few rays, full path)
zpath=np.arange(-0.595,2.0,0.005)
ring.rng=np.random.default_rng(11)
traj={}
for nm,M in (('85',M85),('87',M87)):
    X,V,qm=ring_beam(s,8,M,10e3,sigp=0,sigE=0,w=0)
    az=np.arange(8)*2*np.pi/8
    X[:,0]=0.03*np.cos(az); X[:,1]=0.03*np.sin(az)
    X,V,qm=ring_beam(s,8,M,10e3,sigp=0,sigE=0,w=0)  # regenerate consistent velocities
    rec=track(s,X,V,qm,zpath,0.02/(qm*B0))
    traj[nm]=dict(x=rec[:,:,0].T.round(5).tolist(), y=rec[:,:,1].T.round(5).tolist())
out['traj']=dict(z=zpath.round(4).tolist(), **traj)
# big MC
ring.rng=np.random.default_rng(5)
zp=np.arange(1.2,1.9,0.0025)
mc=run_ring(s,N=20000,zplanes=zp)
r7a=np.hypot(mc['87'][...,0],mc['87'][...,1]); r5a=np.hypot(mc['85'][...,0],mc['85'][...,1])
rms7=np.sqrt(np.nanmean(r7a**2,1)); rms5=np.sqrt(np.nanmean(r5a**2,1)); k=np.nanargmin(rms7)
out['env']=dict(z=zp.round(4).tolist(), rms87=(rms7*1e3).round(3).tolist(), rms85=(rms5*1e3).round(3).tolist(), zbest=float(zp[k]))
r7=r7a[k]; r5=r5a[k]
ra,n7,n5,pur=purity_curve(r7,r5)
out['pur']=dict(a=(ra*1e3).round(3).tolist(), T87=n7.round(5).tolist(), L85=n5.round(6).tolist(), pur=pur.round(6).tolist())
bins=np.linspace(0,25,126)
out['hist']=dict(edges=bins.tolist(), h87=np.histogram(r7*1e3,bins)[0].tolist(), h85=np.histogram(r5*1e3,bins)[0].tolist(), N=20000)
xy={}
for nm in ('85','87'):
    P=mc[nm][k][:2500]; xy[nm]=dict(x=(P[:,0]*1e3).round(3).tolist(), y=(P[:,1]*1e3).round(3).tolist())
out['xy']=xy
# pencil beamlet (single azimuth): Ryan's original picture
ring.rng=np.random.default_rng(9)
pen={}
for nm,M in (('85',M85),('87',M87)):
    X,V,qm=ring_beam(s,1500,M,10e3,naz=1,az0=0.0)
    rec=track(s,X,V,qm,np.array([1.0,zp[k]]),0.02/(qm*B0))
    P=rec[1]; pen[nm]=dict(x=(P[:,0]*1e3).round(3).tolist(), y=(P[:,1]*1e3).round(3).tolist(),
        exit_dir=[float(np.nanmean(rec[0][:,3]/rec[0][:,5])),float(np.nanmean(rec[0][:,4]/rec[0][:,5]))],
        exit_pos=[float(np.nanmean(rec[0][:,0])),float(np.nanmean(rec[0][:,1]))])
out['pencil']=pen
json.dump(out,open('base.json','w'))
print('zbest',zp[k])
for a in [1,1.5,2,2.5,3,4]:
    i=np.argmin(abs(ra-a*1e-3)); print(a,n7[i],n5[i],pur[i])
print(pen['85']['exit_dir'],pen['87']['exit_dir'],pen['85']['exit_pos'],pen['87']['exit_pos'])
