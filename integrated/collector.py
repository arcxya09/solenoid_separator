# Collector end: E2 suppressor ring, knife-edge cone hole (Rb-87), graphite dump (Rb-85), deceleration pocket.
# Axisymmetric electrostatics (sep.Axi) + 3-D tracking in the solenoid fringe field, with the pocket's own space charge.
from sep import *

def geometry(design, a=None, Vland=300.0, gap=None, k=None, style='pocket', Ldrift=None, Rg=None, hgrid=0.2e-3, tout_deg=35.0):
    """electrode masks in collector coordinates (zz = z - zc, metres). k: length scale (1 for 10 kV, 2.5 for 30 kV).
    style 'pocket': cone hole -> 12 mm gap -> cylindrical pocket at V0-Vland (original design).
    style 'sphere': cone hole -> grounded drift tube (Ldrift) -> grounded spherical grid (radius Rg about the focus,
    transparent) -> concentric spherical collector cup at V0-Vland, so the retarding field is radial along the
    (nearly radial) Rb-87 paths."""
    d=DESIGNS[design]; a=a or d['a']; V0=d['E']
    k=k or (1.0 if d['E']<20e3 else 2.5)
    Vp=V0-Vland
    Lc=10e-3*k; tin=0.30; tout=np.tan(np.radians(tout_deg))
    zd=Lc; td=3e-3                    # dump plate front face, thickness
    zE2=d['zE2']-d['zc']
    def cone(R,Z):    # thin-walled knife-edge cone, tip at zz=0, merges into the dump plate
        rin=a+tin*Z; rout=a+0.4e-3+tout*Z
        return (Z>=0)&(Z<=zd+td)&(R>=rin)&(R<=np.maximum(rout,rin+0.4e-3))
    Rd=d.get('col_Rd',36e-3*k); Rh=d.get('col_Rh',38e-3*k)    # dump outer radius, housing radius
    def dump(R,Z):    # graphite dump plate (ground), with the continued hole
        return (Z>=zd)&(Z<=zd+td)&(R>=a+tin*Z)&(R<=Rd)
    rE2=d.get('rE2',20e-3)
    def e2(R,Z):
        return (Z>=zE2)&(Z<=zE2+3e-3)&(R>=rE2)&(R<=rE2+10e-3)
    if style=='pocket':
        g=gap or 12e-3*k
        rpk=max(15e-3*k,2.4*a); zpk=zd+td+g; Lpk=100e-3*k; wpk=3e-3; rl=max(10e-3*k,1.6*a)
        def pocket(R,Z):
            wall=(Z>=zpk)&(Z<=zpk+Lpk)&(R>=rpk)&(R<=rpk+wpk)
            lip=(Z>=zpk)&(Z<=zpk+wpk)&(R>=rl)&(R<=rpk)
            bot=(Z>=zpk+Lpk)&(Z<=zpk+Lpk+wpk)&(R<=rpk+wpk)
            return wall|lip|bot
        zend=zpk+Lpk+wpk
        def housing(R,Z):
            return ((R>=Rh)&(R<=Rh+2e-3*k)&(Z>=zd))|((Z>=zend+9e-3*k)&(Z<=zend+11e-3*k))
        parts=[('cone',cone,0.0),('dump',dump,0.0),('pocket',pocket,Vp),('housing',housing,0.0),('E2',e2,-300.0)]
        P=dict(a=a,Vp=Vp,zd=zd,td=td,zpk=zpk,Lpk=Lpk,rpk=rpk,rl=rl,g=g,zbot=zpk+Lpk)
    else:
        Ld=Ldrift or 30e-3*k; Rg=Rg or (zd+td+Ld)        # grid on the sphere through the end of the drift tube
        gapc=gap or 25e-3*k; Rc=Rg+gapc; th=np.radians(d.get('col_th',17))  # cup half-angle (> 86 mrad x margin, < tube)
        rt=(a+tin*(zd+td))+2e-3*k                       # drift-tube inner radius
        def tube(R,Z): return (Z>=zd+td)&(Z<=Rg*np.cos(np.arcsin(min(rt/Rg,0.99))))&(R>=rt)&(R<=rt+2e-3)
        rho=lambda R,Z: np.hypot(R,Z)
        def grid(R,Z):  # one-cell-thick transparent shell, ground
            return (np.abs(rho(R,Z)-Rg)<=0.6*hgrid)&(Z>0)&(R<=Rg*np.sin(th)+2e-3*k)
        def cup(R,Z):
            sh=(rho(R,Z)>=Rc)&(rho(R,Z)<=Rc+3e-3)&(Z>0)&(R<=(Rc+3e-3)*np.sin(th))
            side=(R>=Rc*np.sin(th))&(R<=Rc*np.sin(th)+3e-3)&(Z>=Rg*np.cos(th))&(Z<=Rc*np.cos(th)+3e-3)
            return sh|side
        def shield(R,Z):  # grounded can closing the grid edge to the dump plate
            return (R>=Rg*np.sin(th)+2e-3*k)&(R<=Rg*np.sin(th)+4e-3*k)&(Z>=zd+td)&(Z<=Rg*np.cos(th))
        zend=Rc+3e-3
        def housing(R,Z):
            return ((R>=max((Rc+6e-3)*np.sin(th)+6e-3*k,Rh))&(R<=max((Rc+6e-3)*np.sin(th)+8e-3*k,Rh+2e-3*k))&(Z>=zd))|((Z>=zend+9e-3*k)&(Z<=zend+11e-3*k))
        parts=[('cone',cone,0.0),('dump',dump,0.0),('pocket',cup,Vp),('housing',housing,0.0),('E2',e2,-300.0),
               ('tube',tube,0.0),('shield',shield,0.0),('grid',grid,0.0)]
        P=dict(a=a,Vp=Vp,zd=zd,td=td,Rg=Rg,Rc=Rc,Ld=Ld,rt=rt,gap=gapc,zbot=Rc*np.cos(th))
    P.update(tout=tout,Rd=Rd,Rh=Rh,style=style,k=k,Vland=Vland,zE2=zE2,zend=zend,Lc=Lc,tin=tin)
    return P,parts

TRANSPARENT=('grid',)

class Collector:
    def __init__(s,design,h=0.2e-3,**gk):
        d=DESIGNS[design]; s.design=design; s.d=d
        s.P,parts=geometry(design,hgrid=h,**gk); P=s.P
        z0=P['zE2']-0.06; z1=P['zend']+14e-3*P['k']; rmax=max(42e-3*P['k'],P.get('Rc',0)*0.5+20e-3*P['k'],P['Rh']+6e-3*P['k'])
        s.zc=d['zc']
        # Axi works in absolute z
        el=[(lambda R,Z,f=f: f(R,Z-s.zc),v) for nm,f,v in parts]
        s.ax=Axi(rmax,z0+s.zc,z1+s.zc,h,h,electrodes=el,Vwall=0.0,zend=(0.0,0.0))
        s.lab=np.zeros(s.ax.fix.shape,int)
        for i,(nm,f,v) in enumerate(parts):
            if nm not in TRANSPARENT: s.lab[f(s.ax.R,s.ax.Z-s.zc)]=i+1
        s.names=['']+[p[0] for p in parts]
        s.rho=None; s.ax.solve()
        s.ax0=None; s._el0=[(f,0.0) for f,v in el]; s._box=(rmax,z0+s.zc,z1+s.zc,h)
    def sc_field(s,rho):
        # space-charge part with all electrodes grounded (superposition with the vacuum solution)
        if s.ax0 is None:
            rmax,z0,z1,h=s._box; s.ax0=Axi(rmax,z0,z1,h,h,electrodes=s._el0,Vwall=0.0,zend=(0.0,0.0))
        s.ax0.solve(rho); return s.ax0.Efield()
    def label(s,r,z):
        i=np.clip(np.round(r/s.ax.hr).astype(int),0,len(s.ax.r)-1); j=np.clip(np.round((z-s.ax.z[0])/s.ax.hz).astype(int),0,len(s.ax.z)-1)
        out=(r>s.ax.r[-1])|(z>s.ax.z[-1])
        return np.where(out,4,s.lab[i,j])
    def track(s,B,X,V,qm,m,I_each=None,Efield=(),dl=0.1e-3,nmax=60000,deposit=False):
        """track until each ion hits an electrode or leaves upstream; returns fate dict.
        I_each: current per macro-particle (A) for space-charge deposition in the pocket region."""
        X=X.copy(); V=V.copy(); N=len(X); alive=np.ones(N,bool)
        hit=np.zeros(N,int); KE=np.full(N,np.nan); Xh=np.full((N,3),np.nan); Vh=np.full((N,3),np.nan)
        v0=np.linalg.norm(V,axis=1).mean(); dt=dl/v0
        Ec=s.ax.Efield(); ax=s.ax
        rho=np.zeros(ax.fix.shape) if deposit else None
        z0=ax.z[0]
        for n in range(nmax):
            a=np.where(alive)[0]
            if len(a)==0: break
            Xa=X[a]; Va=V[a]; q=qm[a][:,None]*dt/2
            bx,by,bz=B.field(Xa[:,0],Xa[:,1],Xa[:,2]); Bv=np.stack([bx,by,bz],1)
            E=Ec(Xa)
            for Ef in Efield: E=E+Ef(Xa)
            Va=Va+q*E; t=q*Bv; s2=2*t/(1+np.sum(t*t,1,keepdims=True)); Va=Va+np.cross(Va+np.cross(Va,t),s2); Va=Va+q*E
            Xn=Xa+Va*dt; rn=np.hypot(Xn[:,0],Xn[:,1])
            if deposit:
                fr=np.clip(rn/ax.hr,0,len(ax.r)-1.001); fz=np.clip((Xn[:,2]-z0)/ax.hz,0,len(ax.z)-1.001)
                i=fr.astype(int); j=fz.astype(int); u=fr-i; w=fz-j; ww=I_each*dt*(Xn[:,2]>s.zc)
                for di,dj,f in ((0,0,(1-u)*(1-w)),(1,0,u*(1-w)),(0,1,(1-u)*w),(1,1,u*w)): np.add.at(rho,(i+di,j+dj),ww*f)
            lb=s.label(rn,Xn[:,2]); up=Xn[:,2]<z0
            stop=(lb>0)|up
            if stop.any():
                ia=a[stop]; hit[ia]=np.where(up[stop],-1,lb[stop]); Xh[ia]=Xn[stop]; Vh[ia]=Va[stop]
                KE[ia]=0.5*m*np.sum(Va[stop]**2,1)/qe; alive[ia]=False
            X[a]=Xn; V[a]=Va
        if deposit:
            vol=2*np.pi*np.maximum(ax.r,ax.hr/4)*ax.hr*ax.hz; vol[0]=np.pi*(ax.hr/2)**2*ax.hz
            rho=rho/vol[:,None]
        return dict(hit=hit,KE=KE,X=Xh,V=Vh,alive=alive),rho

def start_from(rec,zp,zstart):
    k=int(np.argmin(abs(zp-zstart))); P=rec[k]; ok=np.isfinite(P[:,0]); return P[ok,:3].copy(),P[ok,3:].copy()

def run_collector(design,recs,zp,B,col=None,Isep=0.0,fpocket=0.0,iters=4,Efield=(),zstart=None,nmax=60000,gk=None):
    """recs: transport records (both isotopes). Isep: total beam current (A) for the pocket space charge
    (Rb-87 share); fpocket: neutralisation inside the pocket by trapped secondary electrons."""
    col=col or Collector(design,**(gk or {})); d=DESIGNS[design]
    zstart=zstart or col.ax.z[0]+0.002
    res={}
    X7,V7=start_from(recs['87'],zp,zstart); X5,V5=start_from(recs['85'],zp,zstart)
    q7=np.full(len(X7),qe/(M87*amu)); q5=np.full(len(X5),qe/(M85*amu))
    I7=Isep*AB['87']/len(recs['87'][0]) if Isep>0 else 0.0
    Esc=[]
    for it in range(iters if Isep>0 else 1):
        f7,rho=col.track(B,X7,V7,q7,M87*amu,I_each=I7,Efield=list(Efield)+Esc,deposit=Isep>0,nmax=nmax)
        if Isep>0:
            rho=rho*(1-fpocket); col.rho=rho if col.rho is None else 0.5*rho+0.5*col.rho
            Esc=[col.sc_field(col.rho)]
    f5,_=col.track(B,X5,V5,q5,M85*amu,Efield=list(Efield),nmax=nmax)
    return col,f7,f5

def summarize(col,f7,f5,N7tot,N5tot):
    nm=col.names; P=col.P
    def frac(f,lab,N): return float(np.sum(f['hit']==lab)/N)
    out={}
    for iso,f,N in (('87',f7,N7tot),('85',f5,N5tot)):
        o={('upstream' if l==-1 else nm[l]):frac(f,l,N) for l in sorted(set(f['hit'].tolist()))}
        o['lost_before']=float(1-len(f['hit'])/N)
        out[iso]=o
    pk=f7['hit']==nm.index('pocket')
    Z=f7['X'][:,2]-col.zc
    bottom=pk&(Z>P['zbot']-1e-3)
    out['87']['reflected']=float(np.sum(f7['hit']==-1)/N7tot)
    if pk.any():
        KEa=f7['KE'][pk]; Va=f7['V'][pk]
        nrm=np.degrees(np.arctan2(np.hypot(Va[:,0],Va[:,1]),Va[:,2]))
        out['87']['all_KE']=[float(np.min(KEa)),float(np.median(KEa)),float(np.max(KEa))]
        out['87']['all_angle_deg']=[float(np.median(nrm)),float(np.percentile(nrm,99))]
    out['87']['pocket_bottom']=float(bottom.sum()/N7tot)
    if bottom.any():
        KE=f7['KE'][bottom]; V=f7['V'][bottom]; r=np.hypot(f7['X'][bottom,0],f7['X'][bottom,1])
        ang=np.degrees(np.arctan2(np.hypot(V[:,0],V[:,1]),V[:,2]))
        out['87']['land_KE']=[float(np.min(KE)),float(np.median(KE)),float(np.max(KE))]
        out['87']['land_angle_deg']=[float(np.median(ang)),float(np.percentile(ang,99))]
        out['87']['land_r_mm']=[float(np.median(r)*1e3),float(np.percentile(r,99)*1e3)]
    d5=f5['hit']==nm.index('dump'); c5=f5['hit']==nm.index('cone')
    r5=np.hypot(f5['X'][d5|c5,0],f5['X'][d5|c5,1])
    if len(r5): out['85']['land_r_mm']=[float(np.percentile(r5,1)*1e3),float(np.median(r5)*1e3),float(np.percentile(r5,99)*1e3)]
    out['85']['in_pocket']=frac(f5,nm.index('pocket'),N5tot)
    return out
