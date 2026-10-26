# Ring surface-ionization source: one parameterised geometry for both designs.
# Lengths of the extraction optics scale with `s` (ring radius R set separately); voltages scale with V0.
from es import *
import pickle, json, sys, os
LEN=('w','h','d','sl','t','zl','Ll','gl','sll','rnd','dr','rhalf','zlo','zhi','d2','sl2','t2')
REV=dict(w=1e-3,h=1e-3,flare=67.5,d=20e-3,sl=4e-3,t=3e-3,zl=32e-3,Ll=6e-3,gl=10e-3,sll=8e-3,rnd=1.4e-3,
         dr=5e-5,rhalf=10e-3,zlo=-4e-3,zhi=90e-3)          # revised (section 14) electrodes at s=1
# v3 (accel-decel, review B2/C3): suppressor -2.5 kV at 10 mm, ground at 16 mm, positive einzel +5.8 kV at 30 mm
V3=dict(w=1e-3,h=1e-3,flare=67.5,d=10e-3,d2=6e-3,t2=3e-3,sl2=4e-3,sl=6e-3,t=3e-3,zl=30e-3,Ll=6e-3,gl=10e-3,sll=8e-3,
        rnd=1.4e-3,dr=5e-5,rhalf=10e-3,zlo=-4e-3,zhi=90e-3)
def geom_kw(R=0.03,V0=10e3,Vl=-16e3,s=1.0,w=None,v3=False,Vmid=None):
    c={k:(v*s if k in LEN else v) for k,v in (V3 if v3 else REV).items()}
    if w is not None: c['w']=w
    c.update(R=R,V0=V0,Vl=Vl)
    if v3: c['Vmid']=-2500.0*V0/10e3 if Vmid is None else Vmid
    return c
def zout(s): return 0.085*s
def surf_field(g):
    E=np.hypot(g.Er,g.Ez); res={}
    for m,v,nm in g.parts:
        mm=np.zeros_like(m)
        for a,b in ((1,0),(-1,0),(0,1),(0,-1)): mm|=np.roll(m,(a,b),(0,1))
        sel=(~g.fix)&mm; sel[:,:3]=False; sel[:,-3:]=False; sel[:3,:]=False; sel[-3:,:]=False
        res[nm]=float(E[sel].max()/1e6) if sel.any() else 0.0
    return res
def run(I,N=2000,iters=6,seed=2,T=1400.0,save=None,**gk):
    c=geom_kw(**gk); s=gk.get('s',1.0)
    g=Geom(**c); g.solve(); sf=surf_field(g)
    X,V,m=emit(g,N,T=T,rng=np.random.default_rng(seed))
    rho=None
    for it in range(iters if I>0 else 1):
        out,alive,rnew=track(g,X,V,m,zout=zout(s),I=I,deposit=I>0)
        if I>0:
            rho=rnew if rho is None else 0.5*rho+0.5*rnew; g.solve(rho)
    out,alive,_=track(g,X,V,m,zout=zout(s))
    st,_=beam_stats(out,g); st={k:float(v) for k,v in st.items()}
    ok=np.isfinite(out[:,0]); x,y,z,vx,vy,vz=out[ok].T; rr=np.hypot(x,y)
    st.update(I=I,trans=float(ok.mean()),tp_rms=float(np.std((x*vy-y*vx)/rr/vz)),
              E_rms=float(np.std(0.5*m*(vx**2+vy**2+vz**2)/qe)),surf=sf,cfg=c,
              phimin=float(g.interp(g.phi,np.full(len(g.z),g.R),g.z)[0][g.z>c['zl']].min()))
    if save:
        pickle.dump(dict(out=out,cfg=c,I=I,m=m,T=T,stats=st,zout=zout(s),
                         phi=g.phi[::4,::4],r=g.r[::4],z=g.z[::4],solid=g.solid[::4,::4]),open(save,'wb'))
    return st
def job(a):
    a=dict(a); I=a.pop('I'); return run(I,**a)
if __name__=='__main__':
    from multiprocessing import Pool
    cf=json.load(open(sys.argv[1])); np_=int(sys.argv[3]) if len(sys.argv)>3 else 4
    with Pool(np_) as p: res=p.map(job,cf)
    json.dump(res,open(sys.argv[2],'w'),indent=1)
    for r,c in zip(res,cf):
        print({k:v for k,v in c.items() if k!='save'},'tr %.3f rp %.2f mean %.2f tp %.2f wid %.2f eps %.3f Emax %.2f phimin %.0f'%(
            r['trans'],r['rp_rms']*1e3,r['rp_mean']*1e3,r['tp_rms']*1e3,r['rwid_rms']*1e3,r['eps_rms']*1e6,max(r['surf'].values()),r['phimin']))
