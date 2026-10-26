# Integrated start-to-end model of the superconducting-solenoid Rb isotope separator.
# One parameter set per design (DESIGNS), one field model, one tracker, one space-charge solver.
import numpy as np, os, pickle, json
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.special import ellipk, ellipe
mu0=4e-7*np.pi; qe=1.602176634e-19; amu=1.66053907e-27; eps0=8.8541878128e-12; me=9.1093837e-31
M85=84.911789738; M87=86.909180531
AB={'85':0.7217,'87':0.2783}; MASS={'85':M85,'87':M87}
HERE=os.path.dirname(os.path.abspath(__file__))
CACHE=os.environ.get('SEPCACHE',HERE+'/cache')

# ---------------------------------------------------------------- designs
DESIGNS={
 'base': dict(E=10e3, R=0.03, s=1.0, Vl=5.8e3, B0=2.80, I=1e-4,   # v3 accel-decel source, 0.1 mA unneutralised
              coil=dict(L=1.0,a_in=0.10,a_out=0.13,nz=250,nr=4,rmax=0.09,zmin=-1.0,zmax=2.5,dr=1e-3,dz=2e-3),
              zs=-0.6, zc=1.365, a=5.9e-3, zE2=1.28, rpipe=0.075,
              steer=[(-0.575,0.05),(-0.275,0.05)]),
 'hc':   dict(E=30e3, R=0.05, s=2.0, Vl=-48e3, B0=2.85, I=1e-3,   # v2-type source x2, 1 mA unneutralised
              coil=dict(L=1.7,a_in=0.15,a_out=0.18,nz=340,nr=3,rmax=0.13,zmin=-1.0,zmax=3.6,dr=1e-3,dz=2.5e-3),
              zs=-0.6, zc=2.25, a=8.5e-3, zE2=2.15, rpipe=0.11,
              steer=[(-0.575,0.05),(-0.275,0.05)]),
}

import copy
for zs in (-0.3,-0.4):
    v=copy.deepcopy(DESIGNS['base']); v['zs']=zs; v['steer']=[(zs+0.025,0.05),(-0.2,0.05)] if zs<-0.35 else [(zs+0.025,0.05),(-0.15,0.05)]
    DESIGNS['base_zs%g'%(-zs)]=v
v=copy.deepcopy(DESIGNS['base']); v['coil'].update(a_in=0.06,a_out=0.09); DESIGNS['base_bore120']=v
# expert-review working point (R0=60 mm, two axis crossings in the coil), coil C3; 0.5 mA total
v=copy.deepcopy(DESIGNS['base']); v.update(R=0.06,B0=1.9129,I=5e-4,a=9.9e-3,rpipe=0.11)
v['coil']=dict(L=1.0,a_in=0.13,a_out=0.16,nz=250,nr=4,rmax=0.12,zmin=-1.0,zmax=2.5,dr=1e-3,dz=2e-3); DESIGNS['wp60']=v
v=copy.deepcopy(v); v.update(B0=1.894,I=0.0,a=6.9e-3); DESIGNS['wp60_0']=v
# machine B (Ryan's 10 mA report machine): wp60 scaled by sqrt(3) in length, E x3, same B
_k=np.sqrt(3)
DESIGNS['mB']=dict(E=30e3,R=0.104,s=3.0,Vl=17.4e3,B0=1.92,I=6e-3,rE2=0.13,col_Rd=0.185,col_Rh=0.19,col_th=28,
    coil=dict(L=_k,a_in=0.13*_k,a_out=0.16*_k,nz=433,nr=4,rmax=0.20,zmin=-1.8,zmax=3.2,dr=1e-3,dz=3e-3),
    zs=-0.6*_k,zc=1.365*_k,a=22e-3,zE2=1.365*_k-0.15,rpipe=0.20,steer=[(-1.0,0.05),(-0.5,0.05)],sc_h=(0.35e-3,8e-3))

# ---------------------------------------------------------------- solenoid field map
def loop_field(a, z0, r, z):
    zz=z-z0; r=np.maximum(r,1e-9)
    k2=np.minimum(4*a*r/((a+r)**2+zz**2),1-1e-12); K=ellipk(k2); E=ellipe(k2)
    pref=mu0/(2*np.pi*np.sqrt((a+r)**2+zz**2)); den=(a-r)**2+zz**2
    return pref*zz/r*(-K+(a*a+r*r+zz*zz)/den*E), pref*(K+(a*a-r*r-zz*zz)/den*E)

class Solenoid:
    """axisymmetric coil (nz x nr current loops), field map normalised to B0 at the centre"""
    def __init__(s, L=1.0, a_in=0.10, a_out=0.13, nz=250, nr=4, rmax=0.09, zmin=-1.0, zmax=2.5, dr=1e-3, dz=2e-3):
        s.L=L; s.rg=np.arange(0,rmax+dr/2,dr); s.zg=np.arange(zmin,zmax+dz/2,dz); s.dr=dr; s.dz=dz
        R,Z=np.meshgrid(s.rg,s.zg,indexing='ij'); Br=np.zeros_like(R); Bz=np.zeros_like(R)
        for zc in (np.arange(nz)+0.5)/nz*L:
            for ac in a_in+(np.arange(nr)+0.5)/nr*(a_out-a_in):
                br,bz=loop_field(ac,zc,R,Z); Br+=br; Bz+=bz
        c=Bz[0,np.argmin(abs(s.zg-L/2))]; s.Br0=Br/c; s.Bz0=Bz/c
        s.NI_per_T=1/c                                   # total ampere-turns per tesla
        # flux psi(r,z)=int_0^r Bz 2 pi r dr  (for Busch canonical momentum)
        s.psi0=np.concatenate([np.zeros((1,len(s.zg))),np.cumsum(0.5*(s.Bz0[1:]+s.Bz0[:-1])*np.pi*(R[1:]**2-R[:-1]**2),0)],0)
        s.set_B(1.0)
    def set_B(s,B0): s.B0=B0; s.Br=s.Br0*B0; s.Bz=s.Bz0*B0
    def psi(s,r,z):
        fr=np.clip(r/s.dr,0,len(s.rg)-1.001); fz=np.clip((z-s.zg[0])/s.dz,0,len(s.zg)-1.001)
        i=fr.astype(int); j=fz.astype(int); u=fr-i; w=fz-j; P=s.psi0
        return s.B0*(P[i,j]*(1-u)*(1-w)+P[i+1,j]*u*(1-w)+P[i,j+1]*(1-u)*w+P[i+1,j+1]*u*w)
    def field(s,x,y,z):
        r=np.hypot(x,y)
        fr=np.clip(r/s.dr,0,len(s.rg)-1.001); fz=np.clip((z-s.zg[0])/s.dz,0,len(s.zg)-1.001)
        i=fr.astype(int); j=fz.astype(int); u=fr-i; w=fz-j
        def bil(A): return A[i,j]*(1-u)*(1-w)+A[i+1,j]*u*(1-w)+A[i,j+1]*(1-u)*w+A[i+1,j+1]*u*w
        br=bil(s.Br); bz=bil(s.Bz); out=(z<s.zg[0])|(z>s.zg[-1])
        br=np.where(out,0,br); bz=np.where(out,0,bz); rs=np.where(r>0,r,1)
        return br*x/rs, br*y/rs, bz
    def onaxis_dBdz(s): return np.gradient(s.Bz0[0],s.dz)

def get_sol(design, B0=None):
    d=DESIGNS[design] if isinstance(design,str) else design
    os.makedirs(CACHE,exist_ok=True)
    fn=CACHE+'/sol_%s.pkl'%('_'.join('%s%g'%(k,v) for k,v in sorted(d['coil'].items())))
    if os.path.exists(fn): s=pickle.load(open(fn,'rb'))
    else: s=Solenoid(**d['coil']); pickle.dump(s,open(fn,'wb'))
    s.set_B(d['B0'] if B0 is None else B0); return s

class ErrField:
    """solenoid with installation and field-quality errors.
    dx,dy: coil offset (m); tx,ty: tilt about y / about x (rad), around coil centre;
    b1x,b1y: transverse dipole error (fraction of B0) in the coil body;
    b2,a2: normal/skew quadrupole error, fraction of B0 at reference radius r0, in the coil body;
    the body envelope is Bz(z)/B0 so the perturbation fades in the fringe."""
    def __init__(s,sol,dx=0.,dy=0.,tx=0.,ty=0.,b1x=0.,b1y=0.,b2=0.,a2=0.,r0=0.03):
        s.sol=sol; s.zc=sol.L/2; s.__dict__.update(dx=dx,dy=dy,tx=tx,ty=ty,b1x=b1x,b1y=b1y,b2=b2,a2=a2,r0=r0)
        s.zg=sol.zg; s.Bz=sol.Bz; s.B0=sol.B0; s.L=sol.L
    def field(s,x,y,z):
        X=x-s.dx; Y=y-s.dy; Z=z-s.zc
        ct,st=np.cos(s.tx),np.sin(s.tx); cu,su=np.cos(s.ty),np.sin(s.ty)
        x1=ct*X-st*Z; z1=st*X+ct*Z; y2=cu*Y-su*z1; z2=su*Y+cu*z1
        bx,by,bz=s.sol.field(x1,y2,z2+s.zc)
        if s.b1x or s.b1y or s.b2 or s.a2:
            env=np.interp(z2+s.zc,s.sol.zg,s.sol.Bz[0])          # = B0 in the body
            g=env/s.r0
            bx=bx+env*s.b1x+g*(s.b2*y2+s.a2*x1); by=by+env*s.b1y+g*(s.b2*x1-s.a2*y2)
        by1=cu*by+su*bz; bz1=-su*by+cu*bz
        return ct*bx+st*bz1, by1, -st*bx+ct*bz1

def steerers(design,kicks,E=None):
    """kicks: [(thx,thy),(thx,thy)] deflection angles (rad) of ST1, ST2 -> list of uniform field slabs"""
    d=DESIGNS[design]; E=E or d['E']; sl=[]
    for (z0,L),(tx,ty) in zip(d['steer'],kicks): sl.append((z0,L,2*E*tx/L,2*E*ty/L))
    def Ef(X):
        Ex=np.zeros(len(X)); Ey=np.zeros(len(X))
        for z0,L,ex,ey in sl:
            m=(X[:,2]>=z0)&(X[:,2]<z0+L); Ex[m]+=ex; Ey[m]+=ey
        return np.stack([Ex,Ey,0*Ex],1)
    return Ef

# ---------------------------------------------------------------- source phase space -> beam
_src={}
def load_src(fn):
    if fn not in _src:
        d=pickle.load(open(fn,'rb')); o=d['out']; o=o[np.isfinite(o[:,0])]
        x,y,z,vx,vy,vz=o.T; r=np.hypot(x,y)
        # the exit plane still sits in the einzel potential tail: rescale speeds to the energy at ground potential
        KE=0.5*d['m']*(vx**2+vy**2+vz**2)/qe; off=np.median(KE)-d['cfg'].get('V0',10e3)
        f=np.sqrt((KE-off)/KE); vx,vy,vz=vx*f,vy*f,vz*f
        _src[fn]=dict(r=r,vr=(x*vx+y*vy)/r,vt=(x*vy-y*vx)/r,vz=vz,zout=d.get('zout',0.085),I=d['I'],cfg=d['cfg'],m=d['m'])
    return _src[fn]

def beam(sol,design,nm,N,src,seed=0,dx=0.,dy=0.,tx=0.,ty=0.,dR=0.,dE=0.,Bshield=None,src2=None,wmix=0.,halo=0.,hsig=10e-3):
    """ions of isotope nm at the source exit plane zs, from the simulated source phase space (Rb-87).
    Busch: p_theta conserved from the emitter (z_e=zs-zout) where the flux through the ion's start radius is
    psi_e (taken from the field map, i.e. no magnetic shield) or Bshield*pi*r^2 if Bshield is given.
    dx,dy: source offset; tx,ty: source tilt about y/x; dR: ring-radius error; dE: relative extraction-voltage error.
    src2,wmix: linear interpolation of the phase space towards a second source run (same seed) - e.g. einzel error."""
    d=DESIGNS[design]; S=load_src(src); zs=d['zs']; ze=zs-S['zout']
    rng=np.random.default_rng(seed+(87 if nm=='87' else 85)*1000)
    i=rng.integers(0,len(S['r']),N)
    r,vr,vt,vz=[S[k][i].copy() for k in ('r','vr','vt','vz')]
    if src2 is not None and wmix:
        S2=load_src(src2); n=min(len(S2['r']),len(S['r'])); j=np.minimum(i,n-1)
        for a,k in ((r,'r'),(vr,'vr'),(vt,'vt'),(vz,'vz')): a+=wmix*(S2[k][j]-S[k][j])
    re=r+dR                                      # emitter radius ~ exit radius (paraxial ring, < 0.2 mm change)
    r=r+dR
    if halo:   # beam halo: a fraction `halo` of the ions gets an extra random angle (Gaussian, hsig per plane)
        hm=rng.random(N)<halo; nh=hm.sum()
        vr[hm]+=vz[hm]*hsig*rng.standard_normal(nh); vt[hm]+=vz[hm]*hsig*rng.standard_normal(nh)
    M=MASS[nm]; sc=np.sqrt(M87/M*(1+dE)); vr,vt,vz=vr*sc,vt*sc,vz*sc
    qm=qe/(M*amu)
    psi_e=sol.psi(re,np.full(N,ze)) if Bshield is None else Bshield*np.pi*re**2
    psi_s=sol.psi(r,np.full(N,zs))
    vt=vt*1.0+qm*(psi_e-psi_s)/(2*np.pi*r)      # p_theta = m r v_t + q psi/2pi conserved
    az=2*np.pi*rng.random(N)
    ur=np.stack([np.cos(az),np.sin(az)],1); ut=np.stack([-np.sin(az),np.cos(az)],1)
    X=np.concatenate([r[:,None]*ur+np.array([dx,dy]),np.full((N,1),zs)],1)
    V=np.concatenate([vr[:,None]*ur+vt[:,None]*ut,vz[:,None]],1)
    for ang,(a,b) in ((tx,(0,2)),(ty,(1,2))):
        if ang:
            c,s_=np.cos(ang),np.sin(ang); va=c*V[:,a]+s_*V[:,b]; vb=-s_*V[:,a]+c*V[:,b]; V[:,a]=va; V[:,b]=vb
    return X,V,np.full(N,qm)

# ---------------------------------------------------------------- tracker
def track(B, X0, V0, qm, zplanes, dt, Efield=(), rpipe=None, nmax=400000):
    """Boris pusher in B (object with .field) and a list of E-field callables. Records (x,y,z,vx,vy,vz) at the
    first crossing of each z plane; ions with r>rpipe are lost (NaN)."""
    X=X0.copy(); V=V0.copy(); N=len(X); zp=np.asarray(zplanes)
    rec=np.full((len(zp),N,6),np.nan); nxt=np.zeros(N,int); alive=np.ones(N,bool)
    h=(qm*dt/2)[:,None]
    for n in range(nmax):
        a=alive
        bx,by,bz=B.field(X[a,0],X[a,1],X[a,2]); Bv=np.stack([bx,by,bz],1); ha=h[a]; Va=V[a]
        E=None
        for Ef in Efield:
            e=Ef(X[a]); E=e if E is None else E+e
        if E is not None: Va=Va+ha*E
        t=ha*Bv; s2=2*t/(1+np.sum(t*t,1,keepdims=True)); vp=Va+np.cross(Va,t); Va=Va+np.cross(vp,s2)
        if E is not None: Va=Va+ha*E
        Xa=X[a]; Xn=Xa+Va*dt[a][:,None]
        idx=np.where(a)[0]
        while True:                                   # record every plane crossed in this step
            k=nxt[idx]; ok=k<len(zp); zk=np.where(ok,zp[np.minimum(k,len(zp)-1)],np.inf)
            c=ok&(Xn[:,2]>=zk)
            if not c.any(): break
            f=((zk[c]-Xa[c,2])/(Xn[c,2]-Xa[c,2]))[:,None]
            rec[k[c],idx[c],:3]=Xa[c]+f*(Xn[c]-Xa[c]); rec[k[c],idx[c],3:]=Va[c]; nxt[idx[c]]+=1
        X[a]=Xn; V[a]=Va
        if rpipe is not None:
            lost=np.hypot(Xn[:,0],Xn[:,1])>rpipe; alive[idx[lost]]=False
        alive&=nxt<len(zp); alive&=V[:,2]>0
        if not alive.any(): break
    return rec

def ion_dt(design,qm,B0=None,frac=0.02):
    B0=B0 or DESIGNS[design]['B0']; return frac/(qm*B0)

# ---------------------------------------------------------------- axisymmetric Poisson solver
class Axi:
    """phi_rr+phi_r/r+phi_zz=-rho/eps0 on r in [0,rmax] (axis regular), z in [z0,z1].
    electrodes: list of (mask(R,Z)->bool, V). Boundary: phi=Vwall at r=rmax, z ends Dirichlet (zend values) or Neumann."""
    def __init__(s,rmax,z0,z1,hr,hz,electrodes=(),Vwall=0.0,zend=(0.0,0.0)):
        s.r=np.arange(0,rmax+hr/2,hr); s.z=np.arange(z0,z1+hz/2,hz); s.hr=hr; s.hz=hz
        nr,nz=len(s.r),len(s.z); R,Z=np.meshgrid(s.r,s.z,indexing='ij'); s.R=R; s.Z=Z
        fix=np.zeros((nr,nz),bool); val=np.zeros((nr,nz)); s.parts=[]
        for m,v in electrodes:
            mm=m(R,Z); fix|=mm; val[mm]=v; s.parts.append((mm,v))
        s.solid=fix.copy()
        fix[-1,:]=True; val[-1,:]=Vwall
        for j,vv in ((0,zend[0]),(-1,zend[1])):
            if vv is not None: fix[:,j]=True; val[:,j]=np.where(s.solid[:,j],val[:,j],vv)
        s.fix=fix; s.val=val; s.zneu=(zend[0] is None, zend[1] is None)
        idx=np.arange(nr*nz).reshape(nr,nz); I,J=np.meshgrid(np.arange(nr),np.arange(nz),indexing='ij')
        free=~fix; rows=[];cols=[];data=[]
        def add(c,v): rows.append(idx[free]); cols.append(c[free]); data.append(np.broadcast_to(v,I.shape)[free])
        ri=np.where(I>0,s.r[I],1.0)
        cr0=np.where(I>0,-2/hr**2,-4/hr**2)
        add(idx,cr0-2/hz**2)
        add(idx[np.minimum(I+1,nr-1),J],np.where(I>0,1/hr**2+1/(2*ri*hr),4/hr**2))
        add(idx[np.maximum(I-1,0),J],np.where(I>0,1/hr**2-1/(2*ri*hr),0.0))
        jp=np.where(J+1<nz,J+1,J-1); jm=np.where(J-1>=0,J-1,J+1)
        add(idx[I,jp],1/hz**2); add(idx[I,jm],1/hz**2)
        fr=fix.ravel(); rows.append(np.where(fr)[0]); cols.append(np.where(fr)[0]); data.append(np.ones(fr.sum()))
        A=sp.csr_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(nr*nz,nr*nz))
        s.lu=spl.splu(A.tocsc())
    def solve(s,rho=None):
        b=np.zeros(s.fix.shape) if rho is None else -rho/eps0
        b=np.where(s.fix,s.val,b); s.phi=s.lu.solve(b.ravel()).reshape(s.fix.shape)
        Er,Ez=np.gradient(-s.phi,s.hr,s.hz); Er[0]=0; s.Er=Er; s.Ez=Ez; return s.phi
    def _bil(s,A,r,z):
        fr=np.clip(r/s.hr,0,len(s.r)-1.001); fz=np.clip((z-s.z[0])/s.hz,0,len(s.z)-1.001)
        i=fr.astype(int); j=fz.astype(int); u=fr-i; w=fz-j
        return A[i,j]*(1-u)*(1-w)+A[i+1,j]*u*(1-w)+A[i,j+1]*(1-u)*w+A[i+1,j+1]*u*w
    def Efield(s,scale=1.0):
        Er=s.Er*scale; Ez=s.Ez*scale
        def Ef(X):
            r=np.hypot(X[:,0],X[:,1]); z=X[:,2]; out=(z<s.z[0])|(z>s.z[-1])|(r>s.r[-1])
            er=np.where(out,0,s._bil(Er,r,z)); ez=np.where(out,0,s._bil(Ez,r,z)); rs=np.where(r>0,r,1)
            return np.stack([er*X[:,0]/rs,er*X[:,1]/rs,ez],1)
        return Ef

# ---------------------------------------------------------------- space charge
class SpaceCharge:
    """steady-state axisymmetric space charge of the ion beam inside the grounded liner, neutralised by
    fneu(z,r) (fraction of the ion charge cancelled by trapped electrons). Charge deposited from recorded planes."""
    def __init__(s,design,z1=None,hr=0.2e-3,hz=5e-3):
        d=DESIGNS[design]; s.d=d; z1=z1 or d['zc']+0.05
        if 'sc_h' in d: hr,hz=d['sc_h']
        s.ax=Axi(d['rpipe'],d['zs'],z1,hr,hz,zend=(0.0,0.0))
        s.zp=s.ax.z.copy()
    def deposit(s,recs,currents,fneu=None):
        ax=s.ax; rho=np.zeros(ax.fix.shape); hr=ax.hr
        for rec,I in zip(recs,currents):
            N=rec.shape[1]; w=I/N
            r=np.hypot(rec[...,0],rec[...,1]); vz=rec[...,5]
            for j in range(len(ax.z)):
                ok=np.isfinite(r[j])&(vz[j]>0)
                fr=r[j][ok]/hr; i=fr.astype(int); u=fr-i; lam=w/vz[j][ok]
                m=i<len(ax.r)-1; i,u,lam=i[m],u[m],lam[m]
                np.add.at(rho[:,j],i,lam*(1-u)); np.add.at(rho[:,j],i+1,lam*u)
        vol=2*np.pi*np.maximum(ax.r,hr/4)*hr; vol[0]=np.pi*(hr/2)**2   # ring-cell areas (per unit z)
        rho/=vol[:,None]
        if fneu is not None: rho*=1-fneu(ax.R,ax.Z)
        return rho
    def field(s,rho):
        s.ax.solve(rho); return s.ax.Efield()

def fneu_profile(design,f0,z_on=None,z_off=None,ramp=0.02):
    """uniform neutralisation f0 between the einzel barrier (source exit) and the E2 suppressor ring"""
    d=DESIGNS[design]; z_on=d['zs']+0.02 if z_on is None else z_on; z_off=d['zE2'] if z_off is None else z_off
    def f(R,Z): return f0*np.clip((Z-z_on)/ramp,0,1)*np.clip((z_off-Z)/ramp,0,1)
    return f

def run_beam(design,src,N=4000,B=None,B0=None,I=0.0,fneu=None,iters=5,Efield=(),seed=0,zplanes=None,sc=None,relax=0.5,beam_kw=None,verbose=False):
    """start-to-end tracking from the source exit to the collector region, with self-consistent space charge
    of total current I (both isotopes), neutralisation profile fneu. Returns recs dict and the plane z array."""
    d=DESIGNS[design]; sol=get_sol(design,B0); B=B or sol
    beam_kw=beam_kw or {}
    beams={nm:beam(sol,design,nm,N,src,seed=seed,**beam_kw) for nm in ('85','87')}
    if zplanes is None: zplanes=np.arange(d['zc']-0.25,d['zc']+0.2,0.0025)
    rp=d['rpipe']
    if I>0:
        sc=sc or SpaceCharge(design); zall=np.union1d(sc.zp,zplanes); Ef=None; rho=None; conv=[]
        for it in range(iters):
            EL=list(Efield)+([Ef] if Ef else [])
            recs={nm:track(B,*beams[nm],zall,ion_dt(design,beams[nm][2],sol.B0),Efield=EL,rpipe=rp) for nm in beams}
            sel=np.isin(zall,sc.zp)
            rn=sc.deposit([recs['85'][sel],recs['87'][sel]],[I*AB['85'],I*AB['87']],fneu)
            r7=np.hypot(recs['87'][...,0],recs['87'][...,1]); zf=float(zall[np.nanargmin(np.sqrt(np.nanmean(r7**2,1)))])
            conv.append(dict(drho=float(np.abs(rn-rho).sum()/np.abs(rn).sum()) if rho is not None else None,zfocus=zf))
            rho=rn if rho is None else relax*rn+(1-relax)*rho
            Ef=sc.field(rho)
            if verbose: print('SC iter',it,conv[-1:] )
        EL=list(Efield)+[Ef]
        recs={nm:track(B,*beams[nm],zplanes,ion_dt(design,beams[nm][2],sol.B0),Efield=EL,rpipe=rp) for nm in beams}
        return recs,zplanes,dict(sc=sc,rho=rho,Ef=Ef,conv=conv,phimin=float(sc.ax.phi.max()))
    recs={nm:track(B,*beams[nm],zplanes,ion_dt(design,beams[nm][2],sol.B0),Efield=list(Efield),rpipe=rp) for nm in beams}
    return recs,zplanes,{}

# ---------------------------------------------------------------- figures of merit
def pur(n7,n5): return AB['87']*n7/np.maximum(AB['87']*n7+AB['85']*n5,1e-30)
def best(out,zp,P=0.999,recenter=True,amax=0.04):
    """re-optimised collector: plane of minimum Rb-87 rms, hole centred on the Rb-87 centroid, radius maximising
    T87 at purity>=P"""
    a7=out['87']; r7a=np.hypot(a7[...,0],a7[...,1]); rms=np.sqrt(np.nanmean(r7a**2,1)); k=int(np.nanargmin(rms))
    return at_plane(out,zp,k,P,recenter,amax)
def at_plane(out,zp,k,P=0.999,recenter=True,amax=0.04,cx=0.,cy=0.):
    a7=out['87'][k]; a5=out['85'][k]
    if recenter: cx,cy=np.nanmean(a7[:,0]),np.nanmean(a7[:,1])
    r7=np.hypot(a7[:,0]-cx,a7[:,1]-cy); r5=np.hypot(a5[:,0]-cx,a5[:,1]-cy)
    r7=np.where(np.isfinite(r7),r7,np.inf); r5=np.where(np.isfinite(r5),r5,np.inf)
    ra=np.linspace(0,amax,801); n7=np.searchsorted(np.sort(r7),ra)/len(r7); n5=np.searchsorted(np.sort(r5),ra)/len(r5)
    p=pur(n7,n5); ok=p>=P
    T=float(n7[ok].max()) if ok.any() else 0.0; a=float(ra[ok][np.argmax(n7[ok])]) if ok.any() else 0.0
    f=np.isfinite(a7[:,0])
    return dict(T=T,a=a,z=float(zp[k]),cx=float(cx),cy=float(cy),rms87=float(np.sqrt(np.nanmean(r7[f]**2))),
                r87max=float(np.max(r7[f])) if f.any() else np.nan, r85min=float(np.min(r5)))
def fixed(out,zp,zc,a,cx=0.,cy=0.):
    k=int(np.argmin(abs(zp-zc))); a7=out['87'][k]; a5=out['85'][k]
    n7=np.mean(np.hypot(a7[:,0]-cx,a7[:,1]-cy)<a); n5=np.mean(np.hypot(a5[:,0]-cx,a5[:,1]-cy)<a)
    return dict(T=float(n7),L85=float(n5),pur=float(pur(n7,n5)),cx7=float(np.nanmean(a7[:,0])),cy7=float(np.nanmean(a7[:,1])))
def stage(out,zp,a,P=0.999,dz=0.05,cx=None):
    """collector on an x/y/z stage (z within +-dz of the design plane is free): best over planes of the
    recentred hole of FIXED radius a"""
    bestv=None
    for k in range(len(zp)):
        a7=out['87'][k]; a5=out['85'][k]
        c=(np.nanmean(a7[:,0]),np.nanmean(a7[:,1]))
        n7=np.mean(np.hypot(a7[:,0]-c[0],a7[:,1]-c[1])<a); n5=np.mean(np.hypot(a5[:,0]-c[0],a5[:,1]-c[1])<a)
        p=pur(n7,n5); v=n7 if p>=P else 0.0
        if bestv is None or v>bestv['T']: bestv=dict(T=float(v),pur=float(p),z=float(zp[k]),cx=float(c[0]),cy=float(c[1]))
    return bestv
