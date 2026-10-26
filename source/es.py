# Axisymmetric electrostatic solver + ion tracking for a ring surface-ionization source
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
qe=1.602176634e-19; amu=1.66053907e-27; eps0=8.8541878128e-12; kB=8.617333e-5
M87=86.909180531; M85=84.911789738

class Geom:
    def __init__(s, R=0.030, w=1.0e-3, h=0.5e-3, flare=0.0, d=20e-3, sl=4e-3, t=2e-3, V0=10e3,
                 Vmid=None, d2=None, sl2=None, t2=1.5e-3, rc=None, Vl=None, rnd=0.0, zl=None, Ll=6e-3, gl=6e-3, sll=8e-3, loff=0.0,
                 dr=5e-5, rhalf=10e-3, zlo=-4e-3, zhi=70e-3):
        s.__dict__.update(locals()); del s.__dict__['s']
        s.r=np.arange(R-rhalf, R+rhalf+dr/2, dr); s.z=np.arange(zlo, zhi+dr/2, dr)
        Rr,Zz=np.meshgrid(s.r,s.z,indexing='ij'); x=np.abs(Rr-R)
        fix=np.zeros(Rr.shape,bool); val=np.zeros(Rr.shape); s.parts=[]
        # emitter / Pierce block at V0: solid for z<=0 except slot; slot bottom (emitter) at z=-h
        slot_half=w/2+np.maximum(Zz+h,0)*np.tan(np.radians(flare))
        ze=-h-((w/2)**2-np.minimum(x,w/2)**2)/(2*rc) if rc else -h+0*x
        s.ze=(lambda xx: -h-((w/2)**2-np.minimum(np.abs(xx),w/2)**2)/(2*rc)) if rc else (lambda xx: -h+0*xx)
        emi=(Zz<=0)&~((x<slot_half)&(Zz>ze))
        s.parts.append((emi,V0,'emitter'))
        # optional intermediate electrode (e.g. focusing/decel), then ground plate
        if Vmid is not None:
            m=(Zz>=d)&(Zz<=d+t2)&(x>=sl2/2); s.parts.append((m,Vmid,'mid'))
            dg=d+d2
        else: dg=d
        s.dg=dg
        g=(Zz>=dg)&(Zz<=dg+t)&(x>=sl/2); s.parts.append((g,0.0,'ground1'))
        if Vl is not None:   # annular einzel: ground plate (above) - Vl electrode - ground plate
            xl=np.abs(Rr-R-loff)
            m=(Zz>=zl)&(Zz<=zl+Ll)&(xl>=sll/2); s.parts.append((m,Vl,'einzel'))
            m=(Zz>=zl+Ll+gl)&(Zz<=zl+Ll+gl+t)&(xl>=sll/2); s.parts.append((m,0.0,'ground2'))
        if rnd>0:
            from scipy.ndimage import binary_opening
            k=int(round(rnd/dr)); yy,xx=np.mgrid[-k:k+1,-k:k+1]; st=(yy**2+xx**2)<=k*k
            newp=[]
            for m,v,nm in s.parts:
                # pad so that domain edges (plates extending beyond the grid) are not eroded
                mp=np.pad(m,k+1,mode='edge'); mo=binary_opening(mp,structure=st)[k+1:-k-1,k+1:-k-1]
                newp.append((mo,v,nm))
            s.parts=newp
        for m,v,nm in s.parts: fix|=m; val[m]=v
        top=Zz>=zhi-1e-12; fix|=top; val[top]=0.0
        s.fix=fix; s.val=val; s.solid=fix&~top; s.solidf=s.solid.astype(float)
        s._build()
    def _build(s):
        nr,nz=len(s.r),len(s.z); N=nr*nz; idx=np.arange(N).reshape(nr,nz); h=s.dr
        rows=[];cols=[];data=[]
        I,J=np.meshgrid(np.arange(nr),np.arange(nz),indexing='ij')
        free=~s.fix
        ri=s.r[I]
        def add(a,b,v): rows.append(a[free]); cols.append(b[free]); data.append(v[free] if np.ndim(v) else np.full(free.sum(),v))
        c=idx
        add(c,c,np.full(I.shape,-4.0/h**2))
        # r neighbors with Neumann mirror at edges
        ip=np.where(I+1<nr,I+1,I-1); im=np.where(I-1>=0,I-1,I+1)
        jp=np.where(J+1<nz,J+1,J-1); jm=np.where(J-1>=0,J-1,J+1)
        add(c,idx[ip,J],1/h**2+1/(2*ri*h)); add(c,idx[im,J],1/h**2-1/(2*ri*h))
        add(c,idx[I,jp],np.full(I.shape,1/h**2)); add(c,idx[I,jm],np.full(I.shape,1/h**2))
        fr=s.fix.ravel(); 
        rows.append(np.where(fr)[0]); cols.append(np.where(fr)[0]); data.append(np.ones(fr.sum()))
        A=sp.csr_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(N,N))
        s.lu=spl.splu(A.tocsc()); s.N=N
    def solve(s,rho=None):
        b=np.zeros(s.fix.shape)
        if rho is not None: b=-rho/eps0
        b[s.fix]=s.val[s.fix]
        phi=s.lu.solve(b.ravel()).reshape(s.fix.shape)
        s.phi=phi; Er,Ez=np.gradient(-phi,s.dr,s.dr); s.Er=Er; s.Ez=Ez
        return phi
    def interp(s,A,r,z):
        fr=np.clip((r-s.r[0])/s.dr,0,len(s.r)-1.001); fz=np.clip((z-s.z[0])/s.dr,0,len(s.z)-1.001)
        i=fr.astype(int); j=fz.astype(int); u=fr-i; v=fz-j
        return A[i,j]*(1-u)*(1-v)+A[i+1,j]*u*(1-v)+A[i,j+1]*(1-u)*v+A[i+1,j+1]*u*v, (i,j,u,v)

def emit(g, N, M=M87, T=1400.0, rng=None, wem=None):
    rng=rng or np.random.default_rng(0); m=M*amu; kT=kB*T*qe
    we=g.w if wem is None else wem
    r=g.R+we*(rng.random(N)-0.5); z=g.ze(r-g.R)+1.5*g.dr
    sig=np.sqrt(kT/m)
    vz=sig*np.sqrt(-2*np.log(rng.random(N)))           # flux-weighted half-Maxwellian
    vr=sig*rng.standard_normal(N); vt=sig*rng.standard_normal(N)
    return np.stack([r,np.zeros(N),z],1), np.stack([vr,vt,vz],1), m

def track(g, X, V, m, zout, I=0.0, deposit=False, dtmin=2e-11):
    """track in axisymmetric E field; X,V in (r-direction x, azimuthal y, z) local cartesian at azimuth 0"""
    X=X.copy(); V=V.copy(); N=len(X); qm=qe/m
    alive=np.ones(N,bool); out=np.full((N,6),np.nan); done=np.zeros(N,bool)
    rho=np.zeros(g.fix.shape) if deposit else None
    cellvol=2*np.pi*g.r[:,None]*g.dr*g.dr
    for n in range(200000):
        act=alive&~done
        if not act.any(): break
        x,y,z=X[act].T; r=np.hypot(x,y)
        Er,_=g.interp(g.Er,r,z); Ez,(i,j,u,v)=g.interp(g.Ez,r,z)
        vmag=np.linalg.norm(V[act],axis=1)
        dt=np.clip(0.02*g.dr/np.maximum(vmag,1.0),dtmin,5e-9)
        a=qm*np.stack([Er*x/r,Er*y/r,Ez],1)
        Vn=V[act]+a*dt[:,None]; Xn=X[act]+Vn*dt[:,None]
        if deposit:
            w=(I/N)*dt
            for di,dj,ww in ((0,0,(1-u)*(1-v)),(1,0,u*(1-v)),(0,1,(1-u)*v),(1,1,u*v)):
                np.add.at(rho,(i+di,j+dj),w*ww)
        ids=np.where(act)[0]
        # crossing zout
        c=(Xn[:,2]>=zout)
        if c.any():
            f=((zout-X[act][c,2])/(Xn[c,2]-X[act][c,2]))[:,None]
            P=X[act][c]+f*(Xn[c]-X[act][c]); out[ids[c],:3]=P; out[ids[c],3:]=Vn[c]; done[ids[c]]=True
        X[act]=Xn; V[act]=Vn
        # loss: hit electrode
        rr=np.hypot(Xn[:,0],Xn[:,1]); 
        sm,_=g.interp(g.solidf,rr,Xn[:,2])
        hit=(sm>0.75)|(rr<g.r[0])|(rr>g.r[-1])
        alive[ids[hit&~c]]=False
    if deposit: rho=rho/cellvol
    return out, alive, rho

def beam_stats(out,g):
    ok=np.isfinite(out[:,0])
    x,y,z,vx,vy,vz=out[ok].T
    r=np.hypot(x,y); rp=(x*vx+y*vy)/r/vz; tp=(x*vy-y*vx)/r/vz
    dr=r-r.mean(); k=np.polyfit(dr,rp,1)[0]
    res=dict(n=int(ok.sum()), rmean=r.mean(), rwid_rms=dr.std(), rp_mean=rp.mean(), rp_rms=rp.std(), tp_rms=tp.std(),
             slope=k, rp_uncorr=np.std(rp-k*dr), eps_rms=np.sqrt(np.mean(dr**2)*np.mean((rp-rp.mean())**2)-np.mean(dr*(rp-rp.mean()))**2))
    return res, (r,rp,tp)
