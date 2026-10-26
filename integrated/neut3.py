# Louver / deep-groove beam dump for the sputtered-neutral problem (extends neut2.py).
# Collector coords (m): hole plane z=0 (r<a, counted when crossing +z), cone outer r=r0+kc z (0..Lc),
# liner r=Rl (zlo..Lc), E2 ring disk at zE2, upstream exit z=zlo.  Dump region below z=Lc:
#   annulus face z=Lc for r in [Rd,Rl]  (planar case: whole face r in [rc,Rl])
#   louver box r in [rc,Rd], z in [Lc,Lc+D]: inner wall r=rc, outer wall r=Rd, bottom z=Lc+D,
#   thin conical fins r = r_i + m_i (z-Lc) (two-sided).
# All surfaces: sticking s, diffuse (cosine) re-emission about the normal facing the incoming atom.
import numpy as np
from neut2 import cosdir
def _cone(P,D,A,m,z0,zlo,zhi):
    """ray vs cone r=A+m(z-z0), z in [zlo,zhi]; A,m arrays (F,). returns t (N,F) (inf if none), and Q-normals later"""
    c0=A[None,:]+m[None,:]*(P[:,2:3]-z0); c1=m[None,:]*D[:,2:3]
    A2=(D[:,0]**2+D[:,1]**2)[:,None]; B2=(2*(P[:,0]*D[:,0]+P[:,1]*D[:,1]))[:,None]; C2=(P[:,0]**2+P[:,1]**2)[:,None]
    a=A2-c1**2; b=B2-2*c0*c1; c=C2-c0**2; disc=b*b-4*a*c
    T=np.full(c0.shape,np.inf)
    with np.errstate(divide='ignore',invalid='ignore'):
        sq=np.sqrt(np.where(disc>=0,disc,0))
        for sg in (-1,1):
            t=np.where(np.abs(a)>1e-14,(-b+sg*sq)/(2*a),-c/b)
            z=P[:,2:3]+t*D[:,2:3]; ok=(disc>=0)&(t>1e-10)&(z>=zlo)&(z<=zhi)&(c0+c1*t>=0)
            T=np.where(ok&(t<T),t,T)
    return T
class Geo:
    def __init__(s,a,r0,kc,Lc,Rd,Rl,zlo,zE2,rE2,rE2o,D=0.0,fins=None,rc=None):
        s.__dict__.update(a=a,r0=r0,kc=kc,Lc=Lc,Rd=Rd,Rl=Rl,zlo=zlo,zE2=zE2,rE2=rE2,rE2o=rE2o,D=D)
        s.rc=r0+kc*Lc if rc is None else rc
        s.fins=fins          # (r_top array, slope array) or None
    def hits(s,P,D):
        N=len(P); T=np.full(N,np.inf); S=np.full(N,-1)
        def upd(t,sid):
            m=(t>1e-10)&(t<T); T[m]=t[m]; S[m]=sid
        with np.errstate(divide='ignore',invalid='ignore'):
            louv=s.D>0
            disks=[(0.,0,s.a,0),(s.zlo,0,s.Rl,9),(s.zE2,s.rE2,s.rE2o,3),(s.Lc,s.Rd if louv else s.rc,s.Rl,1)]
            if louv: disks.append((s.Lc+s.D,s.rc,s.Rd,5))
            for z0,rlo,rhi,sid in disks:
                t=(z0-P[:,2])/D[:,2]; Q=P+t[:,None]*D; r=np.hypot(Q[:,0],Q[:,1]); upd(np.where((r>=rlo)&(r<=rhi),t,np.inf),sid)
            cyl=[(s.Rl,s.zlo,s.Lc,2)]
            if louv: cyl+= [(s.rc,s.Lc,s.Lc+s.D,6),(s.Rd,s.Lc,s.Lc+s.D,7)]
            for R,zl,zh,sid in cyl:
                t=_cone(P,D,np.array([R]),np.array([0.]),0.,zl,zh)[:,0]; upd(t,sid)
            t=_cone(P,D,np.array([s.r0]),np.array([s.kc]),0.,0.,s.Lc)[:,0]; upd(t,4)
            if louv and s.fins is not None:
                rt,mf=s.fins; TF=_cone(P,D,rt,mf,s.Lc,s.Lc,s.Lc+s.D); j=np.argmin(TF,1); t=TF[np.arange(N),j]
                m=(t<T); T[m]=t[m]; S[m]=100+j[m]
        Q=P+np.where(np.isfinite(T),T,0)[:,None]*D; r=np.maximum(np.hypot(Q[:,0],Q[:,1]),1e-12)
        nrm=np.zeros((N,3)); nrm[:,2]=1.0
        radial=np.isin(S,[2,6,7]); nrm[radial]=np.stack([Q[radial,0]/r[radial],Q[radial,1]/r[radial],0*r[radial]],1)
        for sid,mm in [(4,None)]+([(100+k,None) for k in range(len(s.fins[0]))] if (s.D>0 and s.fins is not None) else []):
            w=S==sid
            if not w.any(): continue
            kk=s.kc if sid==4 else s.fins[1][sid-100]
            v=np.stack([Q[w,0]/r[w],Q[w,1]/r[w],np.full(w.sum(),-kk)],1); nrm[w]=v/np.linalg.norm(v,axis=1)[:,None]
        flip=np.sum(nrm*D,1)>0; nrm[flip]*=-1          # normal faces the incoming ray
        return T,S,nrm,Q
def transport(g,P,D,w,s,nb=60,rng=None):
    """P,D: launch points/directions. returns boolean mask of atoms entering the hole (w unused)"""
    N=len(P); alive=np.ones(N,bool); got=np.zeros(N,bool)
    for b in range(nb):
        idx=np.where(alive)[0]
        if len(idx)==0: break
        T,S,Nn,Q=g.hits(P[idx],D[idx])
        hole=(S==0)&(D[idx][:,2]>0); got[idx[hole]]=True
        lost=(S==9)|(S<0)|hole; stick=(rng.random(len(idx))<s)&~lost; cont=~lost&~stick
        P[idx]=Q; D[idx[cont]]=cosdir(rng,Nn[cont],cont.sum()); P[idx[cont]]+=1e-7*D[idx[cont]]
        alive[idx[~cont]]=False
    return got
def ion_impacts(g,X,V):
    """straight-line Rb-85 ions from the z=Lc plane into the louver box: impact point, surface id, cos(incidence)"""
    D=V/np.linalg.norm(V,axis=1)[:,None]; P=X+1e-7*D
    T,S,Nn,Q=g.hits(P,D); cth=np.abs(np.sum(Nn*D,1))
    return Q,S,Nn,cth
