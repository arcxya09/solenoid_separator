# working-point study: ring radius R0 and number of internal crossings at fixed collector distance
import sep, numpy as np, json, sys
from sep import *
import os
SRC=os.environ.get('SRC','src_v3.pkl')
COILS={'C1':dict(L=1.0,a_in=0.10,a_out=0.13,nz=250,nr=4,rmax=0.09,zmin=-1.0,zmax=2.5,dr=1e-3,dz=2e-3),
       'C3':dict(L=1.0,a_in=0.13,a_out=0.16,nz=250,nr=4,rmax=0.12,zmin=-1.0,zmax=2.5,dr=1e-3,dz=2e-3)}
def mk(name,R0,B0,coil='C3',E=10e3,zc=1.365,a=5.9e-3):
    c=COILS[coil]; d=dict(E=E,R=R0,s=1.0,Vl=-16e3,B0=B0,coil=c,zs=-0.6,zc=zc,a=a,zE2=zc-0.085,
              rpipe=c['a_in']-0.02,steer=[(-0.575,0.05),(-0.275,0.05)])
    DESIGNS[name]=d; return d
_old_dt=sep.ion_dt
sep.ion_dt=lambda design,qm,B0=None,frac=0.02: frac/(qm*2.8)       # same step for every field level
def run(name,N=2000,I=0.0,zplanes=None,seed=0,beam_kw=None,f0=0.0,iters=5):
    d=DESIGNS[name]; Rs=round(float(np.median(load_src(SRC)['r'])),3); bk=dict(dR=d['R']-Rs); bk.update(beam_kw or {})
    if zplanes is None: zplanes=np.arange(0.9,2.2,0.005)
    fneu=fneu_profile(name,f0) if f0 else None
    return run_beam(name,SRC,N=N,I=I,zplanes=zplanes,seed=seed,beam_kw=bk,fneu=fneu,iters=iters)
def focus(recs,zp):
    r7=np.hypot(recs['87'][...,0],recs['87'][...,1]); rms=np.sqrt(np.nanmean(r7**2,1)); k=int(np.nanargmin(rms)); return zp[k],rms[k]
