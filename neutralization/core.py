import numpy as np
from scipy.special import ellipk, ellipe
mu0=4e-7*np.pi; qe=1.602176634e-19; amu=1.66053907e-27; eps0=8.8541878128e-12
M85=84.911789738; M87=86.909180531

def loop_field(a, z0, r, z):
    """single current loop radius a at z0, I=1 A; returns Br,Bz"""
    zz=z-z0; r=np.maximum(r,1e-9)
    k2=4*a*r/((a+r)**2+zz**2); k2=np.minimum(k2,1-1e-12)
    K=ellipk(k2); E=ellipe(k2)
    pref=mu0/(2*np.pi*np.sqrt((a+r)**2+zz**2))
    den=(a-r)**2+zz**2
    Bz=pref*(K+(a*a-r*r-zz*zz)/den*E)
    Br=pref*zz/r*(-K+(a*a+r*r+zz*zz)/den*E)
    return Br,Bz

class Solenoid:
    def __init__(s, B0=3.0, L=1.0, a_in=0.10, a_out=0.13, nz=250, nr=4,
                 rmax=0.09, zmin=-1.0, zmax=3.5, dr=1e-3, dz=2e-3):
        s.L=L; s.rg=np.arange(0,rmax+dr/2,dr); s.zg=np.arange(zmin,zmax+dz/2,dz)
        R,Z=np.meshgrid(s.rg,s.zg,indexing='ij')
        Br=np.zeros_like(R); Bz=np.zeros_like(R)
        for zc in (np.arange(nz)+0.5)/nz*L:
            for ac in a_in+(np.arange(nr)+0.5)/nr*(a_out-a_in):
                br,bz=loop_field(ac,zc,R,Z); Br+=br; Bz+=bz
        c=Bz[0,np.argmin(abs(s.zg-L/2))]
        s.scale0=1/c; s.Br0=Br/c; s.Bz0=Bz/c; s.set_B(B0)
        s.dr=dr; s.dz=dz
    def set_B(s,B0):
        s.B0=B0; s.Br=s.Br0*B0; s.Bz=s.Bz0*B0
    def onaxis(s): return s.zg, s.Bz[0]
    def field(s,x,y,z):
        r=np.hypot(x,y)
        fr=np.clip(r/s.dr,0,len(s.rg)-1.001); fz=np.clip((z-s.zg[0])/s.dz,0,len(s.zg)-1.001)
        i=fr.astype(int); j=fz.astype(int); u=fr-i; w=fz-j
        def bil(A): return (A[i,j]*(1-u)*(1-w)+A[i+1,j]*u*(1-w)+A[i,j+1]*(1-u)*w+A[i+1,j+1]*u*w)
        br=bil(s.Br); bz=bil(s.Bz)
        out=(z<s.zg[0])|(z>s.zg[-1]); br=np.where(out,0,br); bz=np.where(out,0,bz)
        rs=np.where(r>0,r,1)
        return br*x/rs, br*y/rs, bz

def track(sol, X0, V0, qm, zplanes, dt, nmax=200000, Efield=None, keep_path=False, path_every=10):
    """Boris pusher. X0,V0:(N,3). qm:(N,) charge/mass. Records state at each z plane (first crossing)."""
    X=X0.copy(); V=V0.copy(); N=len(X); zp=np.asarray(zplanes)
    rec=np.full((len(zp),N,6),np.nan); done=np.zeros((len(zp),N),bool)
    dtv=dt[:,None] if np.ndim(dt) else dt
    path=[]
    for n in range(nmax):
        bx,by,bz=sol.field(X[:,0],X[:,1],X[:,2])
        B=np.stack([bx,by,bz],1)
        h=(qm[:,None]*dtv/2)
        if Efield is not None:
            E=Efield(X); V=V+h*E
        t=h*B; s2=2*t/(1+np.sum(t*t,1,keepdims=True))
        vp=V+np.cross(V,t); V=V+np.cross(vp,s2)
        if Efield is not None: V=V+h*E
        Xn=X+V*dtv
        for k,zk in enumerate(zp):
            c=(~done[k])&(X[:,2]<zk)&(Xn[:,2]>=zk)
            if c.any():
                f=((zk-X[c,2])/(Xn[c,2]-X[c,2]))[:,None]
                rec[k,c,:3]=X[c]+f*(Xn[c]-X[c]); rec[k,c,3:]=V[c]; done[k,c]=True
        X=Xn
        if keep_path and n%path_every==0: path.append(X.copy())
        if done[-1].all() or (X[:,2]>zp[-1]).all(): break
    return (rec,np.array(path)) if keep_path else rec

def ion(mass_amu, E_eV):
    m=mass_amu*amu; v=np.sqrt(2*E_eV*qe/m); return m,v,qe/m

import os, pickle
def get_sol(B0=3.0, **kw):
    fn='sol_cache_%s.pkl'%abs(hash(tuple(sorted(kw.items()))))
    if os.path.exists(fn): s=pickle.load(open(fn,'rb'))
    else: s=Solenoid(**kw); pickle.dump(s,open(fn,'wb'))
    s.set_B(B0); return s
