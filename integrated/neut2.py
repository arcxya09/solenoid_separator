# Sputtered Rb-85 neutral transport around the collector (generalised beamline/neut.py to collector.geometry):
# atoms start at the simulated Rb-85 impact points (cone outer surface / dump face), cosine emission, diffuse
# re-emission with sticking probability s on every surface; counted when they cross the hole plane z=0, r<a, going +z.
# Surfaces (collector coords, m): hole disk; cone outer r=r0+kc z (0..Lc); dump face z=Lc (r_c(Lc)..Rd);
# liner r=Rl (zlo..Lc); E2 ring disk z=zE2 (rE2..rE2o, both faces); upstream exit z=zlo.
import numpy as np
def cosdir(rng,n,N):
    u=rng.random(N); ph=2*np.pi*rng.random(N); ct=np.sqrt(1-u); st=np.sqrt(u)
    a=np.where(np.abs(n[:,2:3])<0.9,np.array([[0,0,1.]]),np.array([[1.,0,0]]))
    e1=np.cross(n,a); e1/=np.linalg.norm(e1,axis=1)[:,None]; e2=np.cross(n,e1)
    return ct[:,None]*n+st[:,None]*(np.cos(ph)[:,None]*e1+np.sin(ph)[:,None]*e2)
def make(a,r0,kc,Lc,Rd,Rl,zlo,zE2,rE2,rE2o):
    def hits(P,D):
        N=len(P); T=np.full(N,np.inf); S=np.full(N,-1); Nn=np.zeros((N,3))
        def upd(t,ok,sid,nrm):
            m=ok&(t>1e-12)&(t<T); T[m]=t[m]; S[m]=sid; Nn[m]=nrm[m]
        with np.errstate(divide='ignore',invalid='ignore'):
            for z0,rlo,rhi,sid,nz in ((0.,0,a,0,-1.),(Lc,r0+kc*Lc,Rd,1,-1.),(zE2,rE2,rE2o,3,0.),(zlo,0,Rl,9,1.)):
                t=(z0-P[:,2])/D[:,2]; Q=P+t[:,None]*D; r=np.hypot(Q[:,0],Q[:,1]); ok=(r>=rlo)&(r<=rhi)
                nrm=np.zeros((N,3)); nrm[:,2]=(-np.sign(D[:,2]) if nz==0. else nz); upd(t,ok,sid,nrm)
            A2=D[:,0]**2+D[:,1]**2; B2=2*(P[:,0]*D[:,0]+P[:,1]*D[:,1]); C2=P[:,0]**2+P[:,1]**2-Rl**2
            t=(-B2+np.sqrt(B2*B2-4*A2*C2))/(2*A2); Q=P+t[:,None]*D; ok=(Q[:,2]>zlo)&(Q[:,2]<Lc)
            upd(t,ok,2,np.stack([-Q[:,0]/Rl,-Q[:,1]/Rl,0*t],1))
            A=A2-(kc*D[:,2])**2; B=B2-2*kc*D[:,2]*(kc*P[:,2]+r0); C=P[:,0]**2+P[:,1]**2-(kc*P[:,2]+r0)**2
            disc=B*B-4*A*C
            for sg in (-1,1):
                t=(-B+sg*np.sqrt(disc))/(2*A); Q=P+t[:,None]*D; ok=(disc>=0)&(Q[:,2]>=0)&(Q[:,2]<=Lc)
                rq=np.maximum(np.hypot(Q[:,0],Q[:,1]),1e-12)
                nrm=np.stack([Q[:,0]/rq,Q[:,1]/rq,-kc+0*t],1); nrm/=np.linalg.norm(nrm,axis=1)[:,None]; upd(t,ok,4,nrm)
        return T,S,Nn
    return hits
def run(start_r,start_z,geom,s,N=300000,nb=40,seed=1):
    """start_r,start_z: Rb-85 impact points (collector coords). returns fraction of sputtered atoms entering the hole"""
    rng=np.random.default_rng(seed); a,r0,kc,Lc=geom['a'],geom['r0'],geom['kc'],geom['Lc']
    hits=make(**geom)
    i=rng.integers(0,len(start_r),N); r=start_r[i]; ph=2*np.pi*rng.random(N)
    oncone=start_z[i]<Lc-1e-4; z0=np.where(oncone,(r-r0)/kc,Lc); z0=np.clip(z0,0,Lc)
    P=np.stack([r*np.cos(ph),r*np.sin(ph),z0],1)
    n=np.tile([0,0,-1.],(N,1)); nc=np.stack([np.cos(ph),np.sin(ph),np.full(N,-kc)],1); nc/=np.linalg.norm(nc,axis=1)[:,None]
    n[oncone]=nc[oncone]; D=cosdir(rng,n,N); P+=1e-7*D
    alive=np.ones(N,bool); got=0; first=None
    for b in range(nb):
        idx=np.where(alive)[0]; T,S,Nn=hits(P[idx],D[idx])
        hole=(S==0)&(D[idx][:,2]>0); got+=hole.sum()
        if b==0: first=hole.sum()
        lost=(S==9)|(S<0)|hole
        stick=(rng.random(len(idx))<s)&~lost
        P[idx]=P[idx]+np.where(np.isfinite(T),T,0)[:,None]*D[idx]
        cont=~lost&~stick
        D[idx[cont]]=cosdir(rng,Nn[cont],cont.sum()); P[idx[cont]]+=1e-7*D[idx[cont]]
        alive[idx[~cont]]=False
        if not alive.any(): break
    return float(got/N), float(first/N)
