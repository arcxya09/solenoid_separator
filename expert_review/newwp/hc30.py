import wp, sep, sys
from wp import *
from bscan import stage_any
from multiprocessing import Pool
k=np.sqrt(3)
COIL=dict(L=1.0*k,a_in=0.13*k,a_out=0.16*k,nz=433,nr=4,rmax=0.20,zmin=-1.8,zmax=4.3,dr=1e-3,dz=3e-3)
_init=sep.SpaceCharge.__init__
def _sc(s,design,z1=None,hr=0.35e-3,hz=8e-3): _init(s,design,z1,hr,hz)
sep.SpaceCharge.__init__=_sc
ZC=1.365*k
def mk30(nm,B):
    DESIGNS[nm]=dict(E=30e3,R=0.104,s=3.0,Vl=17.4e3,B0=B,coil=COIL,zs=-0.6*k,zc=ZC,a=0.01,zE2=ZC-0.15,rpipe=0.20,steer=[(-1.0,0.05),(-0.5,0.05)])
sep.ion_dt=lambda design,qm,B0=None,frac=0.02: frac/(qm*2.8)
def job(a):
    B,I,src,N=a; wp.SRC=src; nm='h_%g_%g'%(B,I); mk30(nm,B)
    zp=np.arange(ZC-0.2,ZC+0.2,0.005) if I==0 and N<1000 else np.arange(ZC-0.09,ZC+0.09,0.005)
    recs,zp,_=run_beam(nm,src,N=N,I=I*1e-3,zplanes=zp,seed=1,beam_kw=dict(dR=0.0),iters=4)
    if N<1000: zf,rms=focus(recs,zp); return dict(B=B,I=I,zf=float(zf),rms=float(rms))
    s=stage_any(recs,zp,zc=ZC,dz=0.087); s.update(B=B,I=I,src=src); return s
if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='refine':
        mk30('init',1.9); get_sol('init')
        J=[(round(B,3),0,'src104_30k_3mA.pkl',400) for B in np.arange(1.86,1.94,0.01)]
    else:
        B0=float(sys.argv[2])
        J=[(round(B0*(1+f),4),I,s,1500) for I,s,fs in ((0,'src104_30k_3mA.pkl',(0,)),(3,'src104_30k_3mA.pkl',(0.01,0.02,0.03)),(4,'src104_30k_5mA.pkl',(0.015,0.025,0.035)),(5,'src104_30k_5mA.pkl',(0.02,0.03,0.04))) for f in fs]
    with Pool(4) as p:
        for r in p.imap(job,J): print(r,flush=True)
