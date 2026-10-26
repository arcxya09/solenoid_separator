# Self-consistent space-charge neutralization model for design A (30 keV, R=50 mm, 1.7 m x 2.8 T)
# electrons: magnetized, Boltzmann along each flux tube, n_e = n_ref(psi) exp((phi-Phi*(psi))/Te)
# ions: tracked (Boris) in B + E=-grad(phi); 2D axisymmetric Poisson (r,z), wall r=RW and ends grounded
import numpy as np, pickle, os, sys, json
import scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.ndimage import gaussian_filter
from core import *
from scan import metrics
KW=dict(L=1.7,a_in=0.15,a_out=0.18,nz=340,nr=3,rmax=0.13,zmin=-1.0,zmax=5.0,dr=1e-3,dz=2.5e-3)
E0=30e3; R=0.05; ZS=-0.6
B0=float(os.environ.get('B0','2.8'))
s=get_sol(B0,**KW)
RW=0.13; DR=2.5e-4; DZ=5e-3; ZC=float(os.environ.get('ZC','2.6'))
rg=np.arange(0,RW+DR/2,DR); zg=np.arange(ZS,ZC+DZ/2,DZ); NR,NZ=len(rg),len(zg)
def load_src(fn='src_0p1mA.pkl'):
    d=pickle.load(open(fn,'rb')); o=d['out']; o=o[np.isfinite(o[:,0])]
    x,y,z,vx,vy,vz=o.T; r=np.hypot(x,y); return r,(x*vx+y*vy)/r,(x*vy-y*vx)/r,vz
SRC=load_src()
def beam(M,N,rng,scale_ang=1/np.sqrt(3)):
    r,vr,vt,vz=SRC; i=rng.integers(0,len(r),N)
    r=r[i]-r[i].mean()+R; vr=vr[i]; vt=vt[i]; vz=vz[i]
    v87=np.sqrt(2*E0*qe/(M87*amu)); ar=vr/vz*scale_ang; at=vt/vz*scale_ang
    sc=v87*np.sqrt(M87/M); vz=sc/np.sqrt(1+ar**2+at**2); vr=ar*vz; vt=at*vz
    qm=qe/(M*amu); az=2*np.pi*rng.random(N)
    bz0=np.interp(ZS,s.zg,s.Bz[0]); vt=vt-qm*bz0*r/2
    ur=np.stack([np.cos(az),np.sin(az)],1); ut=np.stack([-np.sin(az),np.cos(az)],1)
    X=np.concatenate([r[:,None]*ur,np.full((N,1),ZS)],1)
    V=np.concatenate([vr[:,None]*ur+vt[:,None]*ut,vz[:,None]],1)
    return X,V,np.full(N,qm)
# flux-tube label: equivalent radius in the uniform field, r_eq=sqrt(2 psi/B0)
def _psi():
    Rf,Zf=np.meshgrid(s.rg,s.zg,indexing='ij')
    ps=np.cumsum(s.Bz*Rf,0)*s.dr  # approx int Bz r dr
    from scipy.interpolate import RegularGridInterpolator
    f=RegularGridInterpolator((s.rg,s.zg),ps)
    Rg,Zg=np.meshgrid(rg,zg,indexing='ij')
    return np.sqrt(2*np.maximum(f((Rg,Zg)),0)/B0)
REQ=_psi()
# Poisson operator
def build_L():
    n=NR*NZ; idx=lambda i,j:i*NZ+j
    rows=[];cols=[];vals=[]
    rp=rg+DR/2; rm=np.maximum(rg-DR/2,0)
    I,J=np.meshgrid(np.arange(NR),np.arange(NZ),indexing='ij')
    inner=(I<NR-1)&(J>0)&(J<NZ-1)
    def add(mask,di,dj,v):
        ii=I[mask];jj=J[mask]; rows.append(idx(ii,jj)); cols.append(idx(ii+di,jj+dj)); vals.append(v[mask] if np.ndim(v) else np.full(mask.sum(),v))
    Rr=rg[I]; m0=inner&(I==0); m1=inner&(I>0)
    add(m0,0,0,-4/DR**2); add(m0,1,0,4/DR**2)
    add(m1,0,0,-(rp[I]+rm[I])/(Rr*DR**2+1e-30)); add(m1,1,0,rp[I]/(Rr*DR**2+1e-30)); add(m1,-1,0,rm[I]/(Rr*DR**2+1e-30))
    add(inner,0,1,1/DZ**2); add(inner,0,-1,1/DZ**2); add(inner,0,0,-2/DZ**2)
    bd=~inner; rows.append(idx(I[bd],J[bd])); cols.append(idx(I[bd],J[bd])); vals.append(np.ones(bd.sum()))
    L=sp.csr_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n))
    return L,inner.ravel()
L,INNER=build_L()
def deposit(recs,ws,smooth=1.0):
    """ion density n_i(r,z) [m^-3] from records at planes zg"""
    n=np.zeros((NR,NZ))
    for rec,w in zip(recs,ws):
        r=np.hypot(rec[...,0],rec[...,1]); vz=rec[...,5]
        for j in range(NZ):
            ok=np.isfinite(r[j])&(vz[j]>0); fr=r[j][ok]/DR; i=np.floor(fr).astype(int); u=fr-i
            m=i<NR-1; np.add.at(n[:,j],i[m],w/vz[j][ok][m]*(1-u[m])); np.add.at(n[:,j],i[m]+1,w/vz[j][ok][m]*u[m])
    area=2*np.pi*np.maximum(rg,DR/4)*DR; area[0]=np.pi*(DR/2)**2
    n=n/area[:,None]/qe
    return gaussian_filter(n,(smooth,smooth*0.5))
BINS=np.arange(0,RW+DR,DR); TB=np.clip(np.digitize(REQ,BINS)-1,0,len(BINS)-2)
def newton(ni,NREF,PH,Te,phi,tol=1e-3,verbose=False):
    Q=qe/eps0; ni_=ni.ravel(); NR_=NREF.ravel(); PH_=PH.ravel()
    for it in range(300):
        ne=NR_*np.exp(np.minimum((phi-PH_)/Te,30))
        F=L@phi+Q*(ni_-ne)*INNER
        d=sla.spsolve((L-sp.diags(Q*ne/Te*INNER)).tocsc(),-F)
        mx=np.abs(d).max(); lim=3*Te
        if mx>lim: d*=lim/mx
        phi=phi+d
        if mx<tol*Te: break
    return phi,it
LNREF=None
ERR=None
def solve(ni,Te,Phis,phi0=None,outer=30,verbose=False):
    """Te [V]; Phis = potential maximum on each beam-threaded tube [V] (scalar or fn of r_eq).
    n_ref(tube) is adjusted until max_z phi on the tube equals Phis (barrier-limited trapping)."""
    PH=Phis(REQ) if callable(Phis) else np.full_like(REQ,float(Phis))
    nb=len(BINS)-1
    nimax=np.zeros(nb); np.maximum.at(nimax,TB.ravel(),ni.ravel())
    has=nimax>1e-7*nimax.max()
    lnref=np.where(has,np.log(np.maximum(nimax,1)),-np.inf)
    global LNREF
    if LNREF is not None: lnref=np.where(has,np.where(np.isfinite(LNREF),LNREF,lnref),-np.inf)
    phi=np.zeros(NR*NZ) if phi0 is None else phi0.ravel().copy()
    for o in range(outer):
        NREF=np.exp(lnref)[TB]
        phi,it=newton(ni,NREF,PH,Te,phi)
        pm=np.full(nb,-np.inf); np.maximum.at(pm,TB.ravel(),(phi-PH.ravel()))
        err=np.where(has,pm,0); 
        if verbose: print('outer',o,'newton',it,'max|phimax-Phi*| %.2f'%np.abs(err[has]).max())
        lnref=np.where(has,lnref+np.clip(err,-6*Te,6*Te)/Te,-np.inf)
        if np.abs(err[has]).max()<0.1*Te: break
    LNREF=lnref; global ERR; ERR=float(np.abs(err[has]).max())
    NREF=np.exp(lnref)[TB]; ne=NREF*np.exp(np.minimum((phi.reshape(NR,NZ)-PH)/Te,30))
    return phi.reshape(NR,NZ),ne,o
def efield(phi):
    Er=-np.gradient(phi,DR,axis=0); Ez=-np.gradient(phi,DZ,axis=1); Er[0]=0
    def Ef(X):
        r=np.hypot(X[:,0],X[:,1]); z=X[:,2]
        fr=np.clip(r/DR,0,NR-1.001); fz=np.clip((z-zg[0])/DZ,0,NZ-1.001)
        i=fr.astype(int); j=fz.astype(int); u=fr-i; w=fz-j
        bil=lambda A:A[i,j]*(1-u)*(1-w)+A[i+1,j]*u*(1-w)+A[i,j+1]*(1-u)*w+A[i+1,j+1]*u*w
        er=bil(Er); ez=bil(Ez); out=(z<zg[0])|(z>zg[-1])|(r>RW); er=np.where(out,0,er); ez=np.where(out,0,ez)
        rs=np.where(r>0,r,1); return np.stack([er*X[:,0]/rs,er*X[:,1]/rs,ez],1)
    return Ef
AB={'85':0.7217,'87':0.2783}
def beams(N=3000):
    out=[]
    for k,(nm,M) in enumerate((('85',M85),('87',M87))):
        X,V,qm=beam(M,N,np.random.default_rng(30+k)); out.append((nm,X,V,qm))
    return out
def run_tracks(bm,Ef):
    return [track(s,X,V,qm,zg,0.02/(qm*B0),Efield=Ef) for nm,X,V,qm in bm]
def evaluate(recs,zmin=1.75):
    sel=zg>zmin
    r7a=np.hypot(recs[1][sel][...,0],recs[1][sel][...,1]); r5a=np.hypot(recs[0][sel][...,0],recs[0][sel][...,1])
    rms7=np.sqrt(np.nanmean(r7a**2,1)); k=np.nanargmin(rms7)
    m=metrics(r7a[k],r5a[k]); m.update(z=float(zg[sel][k]),rms87=float(rms7[k])); return m
