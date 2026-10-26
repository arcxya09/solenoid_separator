# Beamline from ring-source exit to collector: misalignment, steering, space charge, scattering, collector
from core import *
import pickle, os
SRC=os.environ.get('SRCDIR','/mnt/project-files/solenoid-separator/beamline/src/')
AB={'85':0.7217,'87':0.2783}; MASS={'85':M85,'87':M87}
_cache={}
def load_src(fn='src_0p1mA.pkl'):
    if fn not in _cache:
        d=pickle.load(open(SRC+fn,'rb')); o=d['out']; o=o[np.isfinite(o[:,0])]
        x,y,z,vx,vy,vz=o.T; r=np.hypot(x,y)
        _cache[fn]=(r,(x*vx+y*vy)/r,(x*vy-y*vx)/r,vz)
    return _cache[fn]

class Mis:
    """coil displaced by (dx,dy) and tilted by (tx about y, ty about x) around its centre z=L/2"""
    def __init__(s,sol,dx=0.,dy=0.,tx=0.,ty=0.):
        s.sol=sol; s.dx=dx; s.dy=dy; s.tx=tx; s.ty=ty; s.zc=sol.L/2
        s.zg=sol.zg; s.Bz=sol.Bz; s.B0=sol.B0
    def field(s,x,y,z):
        X=x-s.dx; Y=y-s.dy; Z=z-s.zc
        ct,st=np.cos(s.tx),np.sin(s.tx); cu,su=np.cos(s.ty),np.sin(s.ty)
        x1=ct*X-st*Z; z1=st*X+ct*Z          # undo tilt about y
        y2=cu*Y-su*z1; z2=su*Y+cu*z1        # undo tilt about x
        bx,by,bz=s.sol.field(x1,y2,z2+s.zc)
        by1=cu*by+su*bz; bz1=-su*by+cu*bz
        return ct*bx+st*bz1, by1, -st*bx+ct*bz1

def beam(sol,nm,N,src='src_0p1mA.pkl',zs=-0.6,dx=0.,dy=0.,tilt=0.,tilty=0.,Bres=0.,dR=0.,dE=0.,seed=0):
    """ions at the source exit plane zs. Bres: residual field at the emitter (flux through ring = pi r^2 Bres).
    dE: relative energy error (HV error)."""
    rng=np.random.default_rng(seed+(87 if nm=='87' else 85)*1000)
    r,vr,vt,vz=load_src(src); i=rng.integers(0,len(r),N)
    r,vr,vt,vz=r[i]+dR,vr[i],vt[i],vz[i]
    M=MASS[nm]; sc=np.sqrt(M87/M*(1+dE)); vr,vt,vz=vr*sc,vt*sc,vz*sc
    qm=qe/(M*amu)
    az=2*np.pi*rng.random(N)
    bz0=np.interp(zs,sol.zg,sol.Bz[0]); vt=vt+qm*(Bres-bz0)*r/2   # Busch: p_theta = q*psi_emitter/2pi
    ur=np.stack([np.cos(az),np.sin(az)],1); ut=np.stack([-np.sin(az),np.cos(az)],1)
    X=np.concatenate([r[:,None]*ur+np.array([dx,dy]),np.full((N,1),zs)],1)
    V=np.concatenate([vr[:,None]*ur+vt[:,None]*ut,vz[:,None]],1)
    for ang,(a,b) in ((tilt,(0,2)),(tilty,(1,2))):
        if ang:
            c,s_=np.cos(ang),np.sin(ang); va=c*V[:,a]+s_*V[:,b]; vb=-s_*V[:,a]+c*V[:,b]; V[:,a]=va; V[:,b]=vb
    return X,V,np.full(N,qm)

def steer_field(st):
    """st: list of (z0,L,Ex,Ey) uniform transverse field slabs (thin parallel-plate steerers)"""
    def Ef(X):
        E=np.zeros_like(X)
        for z0,L,ex,ey in st:
            m=(X[:,2]>=z0)&(X[:,2]<z0+L); E[m,0]+=ex; E[m,1]+=ey
        return E
    return Ef
def kick_to_E(theta,L,E_eV=10e3):
    """transverse field (V/m) over length L giving angle theta for an ion of kinetic energy E"""
    return 2*E_eV*theta/L

ZP=np.arange(1.0,1.9,0.0025)
def run(sol,N=6000,zp=ZP,Ef=None,**kw):
    out={}
    for nm in ('85','87'):
        X,V,qm=beam(sol,nm,N,**kw)
        out[nm]=track(sol,X,V,qm,zp,0.02/(qm*sol.B0),Efield=Ef)
    return out

def pur(n7,n5): return AB['87']*n7/np.maximum(AB['87']*n7+AB['85']*n5,1e-30)
def best(out,zp=ZP,recenter=True,P=0.999):
    """re-optimised collector: best plane, hole centred on Rb-87 centroid, best radius for purity>=P"""
    r7a=np.hypot(out['87'][...,0],out['87'][...,1]); rms=np.sqrt(np.nanmean(r7a**2,1)); k=np.nanargmin(rms)
    a7=out['87'][k]; a5=out['85'][k]
    cx,cy=(np.nanmean(a7[:,0]),np.nanmean(a7[:,1])) if recenter else (0.,0.)
    r7=np.hypot(a7[:,0]-cx,a7[:,1]-cy); r5=np.hypot(a5[:,0]-cx,a5[:,1]-cy)
    ra=np.linspace(0,0.015,601); n7=np.array([(r7<a).mean() for a in ra]); n5=np.array([(r5<a).mean() for a in ra])
    p=pur(n7,n5); ok=p>=P
    T=float(n7[ok].max()) if ok.any() else 0.0; a=float(ra[ok][np.argmax(n7[ok])]) if ok.any() else 0.0
    return dict(T=T,a=a,z=float(zp[k]),cx=float(cx),cy=float(cy),rms87=float(np.sqrt(np.mean(r7**2))))
def fixed(out,zc,a,zp=ZP,cx=0.,cy=0.):
    """fixed collector: plane zc, hole radius a centred at (cx,cy)"""
    k=np.argmin(abs(zp-zc)); a7=out['87'][k]; a5=out['85'][k]
    n7=np.nanmean(np.hypot(a7[:,0]-cx,a7[:,1]-cy)<a); n5=np.nanmean(np.hypot(a5[:,0]-cx,a5[:,1]-cy)<a)
    c7=(np.nanmean(a7[:,0]),np.nanmean(a7[:,1]))
    return dict(T=float(n7),L85=float(n5),pur=float(pur(n7,n5)),cx7=float(c7[0]),cy7=float(c7[1]))
