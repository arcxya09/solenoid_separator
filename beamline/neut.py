# neutral (sputtered Rb-85) transport MC around the collector: fraction reaching the Rb-87 aperture
# geometry (mm, z downstream, collector tip z=0): aperture disk r<5.9 @z=0 (counts);
# knife-edge cone r=5.9@z=0 -> r=9@z=10 (outer surface); dump face plane z=6, r 9..75 (facing upstream);
# liner cylinder r=75, z -600..6; E2 ring disk z=-85, r 20..40 (both faces); upstream exit z=-600 (lost).
import numpy as np, json, sys
rng=np.random.default_rng(1)
def cosdir(n,N):
    u=rng.random(N); ph=2*np.pi*rng.random(N); ct=np.sqrt(u); st=np.sqrt(1-u)
    # build basis around normals n (N,3)
    a=np.where(np.abs(n[:,2:3])<0.9,np.array([[0,0,1.]]),np.array([[1.,0,0]]))
    e1=np.cross(n,a); e1/=np.linalg.norm(e1,axis=1)[:,None]; e2=np.cross(n,e1)
    return ct[:,None]*n+st[:,None]*(np.cos(ph)[:,None]*e1+np.sin(ph)[:,None]*e2)
def hits(P,D):
    """return (t, surface id, normal) of first intersection"""
    N=len(P); T=np.full(N,np.inf); S=np.full(N,-1); Nn=np.zeros((N,3))
    def upd(t,ok,sid,nrm):
        m=ok&(t>1e-9)&(t<T); T[m]=t[m]; S[m]=sid; Nn[m]=nrm[m] if nrm.ndim==2 else nrm
    with np.errstate(divide='ignore',invalid='ignore'):
        # planes z=z0
        for z0,rlo,rhi,sid,nz in ((0.,0,5.9,0,-1.),(6.,9.,75.,1,-1.),(-85.,20.,40.,3,0.),(-600.,0,75.,9,1.)):
            t=(z0-P[:,2])/D[:,2]; Q=P+t[:,None]*D; r=np.hypot(Q[:,0],Q[:,1]); ok=(r>=rlo)&(r<=rhi)
            if nz==0.: nrm=np.zeros((N,3)); nrm[:,2]=-np.sign(D[:,2])
            else: nrm=np.tile([0,0,nz],(N,1)).astype(float)
            upd(t,ok,sid,nrm)
        # cylinder r=75
        a=D[:,0]**2+D[:,1]**2; b=2*(P[:,0]*D[:,0]+P[:,1]*D[:,1]); c=P[:,0]**2+P[:,1]**2-75.**2
        t=(-b+np.sqrt(b*b-4*a*c))/(2*a); Q=P+t[:,None]*D; ok=(Q[:,2]>-600)&(Q[:,2]<6)
        nrm=np.stack([-Q[:,0]/75,-Q[:,1]/75,0*t],1); upd(t,ok,2,nrm)
        # cone outer surface: r = 5.9 + 0.31 z, z in [0,10]  -> r - k z - r0 = 0
        k=0.31; r0=5.9
        A=a-(k*D[:,2])**2; B=b-2*k*D[:,2]*(k*P[:,2]+r0); C=P[:,0]**2+P[:,1]**2-(k*P[:,2]+r0)**2
        disc=B*B-4*A*C
        for sg in (-1,1):
            t=(-B+sg*np.sqrt(disc))/(2*A); Q=P+t[:,None]*D; ok=(disc>=0)&(Q[:,2]>=0)&(Q[:,2]<=10)
            rq=np.hypot(Q[:,0],Q[:,1]); nrm=np.stack([Q[:,0]/np.maximum(rq,1e-9),Q[:,1]/np.maximum(rq,1e-9),-k+0*t],1); nrm/=np.linalg.norm(nrm,axis=1)[:,None]
            upd(t,ok,4,nrm)
    return T,S,Nn
def run(s,N=400000,nb=30):
    r=np.sqrt(rng.uniform(7.4**2,16.3**2,N)); ph=2*np.pi*rng.random(N)
    oncone=r<9.0; z0=np.where(oncone,(r-5.9)/0.31,6.0)
    P=np.stack([r*np.cos(ph),r*np.sin(ph),z0],1)
    n=np.tile([0,0,-1.],(N,1)); nc=np.stack([np.cos(ph),np.sin(ph),np.full(N,-0.31)],1); nc/=np.linalg.norm(nc,axis=1)[:,None]
    n[oncone]=nc[oncone]; D=cosdir(n,N); P+=1e-6*D
    alive=np.ones(N,bool); got=0; first=None
    for b in range(nb):
        T,S,Nn=hits(P[alive],D[alive]); idx=np.where(alive)[0]
        hole=(S==0)&(D[alive][:,2]>0); got+=hole.sum()
        if b==0: first=hole.sum()
        lost=(S==9)|(S<0)|hole
        stick=(rng.random(len(idx))<s)&~lost
        P[idx]=P[idx]+np.where(np.isfinite(T),T,0)[:,None]*D[idx]
        cont=~lost&~stick
        D[idx[cont]]=cosdir(Nn[cont],cont.sum()); P[idx[cont]]+=1e-6*D[idx[cont]]
        alive[idx[~cont]]=False
        if not alive.any(): break
    return got/N, first/N
out={}
print('fraction landing on cone', np.mean(np.sqrt(rng.uniform(7.4**2,16.3**2,100000))<9))
for s in (0.3,0.5,0.9,0.99,1.0):
    f,f1=run(s); out[str(s)]=f; print('sticking',s,'fraction into aperture %.2e (direct %.1e)'%(f,f1))
json.dump(out,open('neut.json','w'))
