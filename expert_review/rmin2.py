from sep import *
sol=get_sol('base'); zp=np.arange(-0.6,1.5,0.002); d=DESIGNS['base']
N=200; M=M87; qm=qe/(M*amu); v=np.sqrt(2*10e3*qe/(M*amu)); az=np.linspace(0,2*np.pi,N,endpoint=False)
for Bs in (None,0.0):
    R=0.03; ze=d['zs']-0.085
    psi_e=sol.psi(np.array([R]),np.array([ze]))[0] if Bs is None else 0.0
    psi_s=sol.psi(np.array([R]),np.array([d['zs']]))[0]
    vt=qm*(psi_e-psi_s)/(2*np.pi*R); Bse=psi_e/(np.pi*R*R)
    X=np.stack([R*np.cos(az),R*np.sin(az),np.full(N,d['zs'])],1)
    V=np.stack([-vt*np.sin(az),vt*np.cos(az),np.full(N,np.sqrt(v*v-vt*vt))],1)
    rec=track(sol,X,V,np.full(N,qm),zp,ion_dt('base',np.full(N,qm)))
    rr=np.hypot(rec[...,0],rec[...,1]).mean(1)
    idx=[i for i in range(1,len(zp)-1) if rr[i]<rr[i-1] and rr[i]<rr[i+1]]
    vr=(rec[...,0]*rec[...,3]+rec[...,1]*rec[...,4])/np.maximum(np.hypot(rec[...,0],rec[...,1]),1e-12)/rec[...,5]
    print('B_emit=%.4f T'%Bse, 'minima (z, r mm):',[(round(zp[i],3),round(rr[i]*1e3,3)) for i in idx])
    k=np.argmin(abs(zp-1.2)); print('  convergence angle at z=1.2: %.1f mrad'%(abs(vr[k].mean())*1e3), 'analytic coil r_min=%.3f mm'%(Bse*R/2.785*1e3))
