# can an electrostatic quadrupole on the 8-electrode steerer ST2 compensate the coil quadrupole error b2/a2?
from hc30 import *
import sys
RQ=0.14; ZQ=(-0.476,0.087)
def quadE(Vq,th):
    c,s=np.cos(2*th),np.sin(2*th)
    def Ef(X):
        x,y,z=X[:,0],X[:,1],X[:,2]; m=(z>=ZQ[0])&(z<ZQ[0]+ZQ[1])
        # phi = Vq/RQ^2 * [ c(x^2-y^2) + s(2xy) ]
        ex=-Vq/RQ**2*(2*c*x+2*s*y); ey=-Vq/RQ**2*(-2*c*y+2*s*x)
        return np.stack([np.where(m,ex,0),np.where(m,ey,0),0*x],1)
    return Ef
def job(a):
    I,src,B,b2,a2,Vq,th,N=a; wp.SRC=src; nm='q_%g'%B; mk30(nm,B); sol=get_sol(nm)
    Bf=ErrField(sol,r0=0.104,b2=b2,a2=a2) if (b2 or a2) else sol
    zp=np.arange(2.2,2.55,0.005)
    EL=[quadE(Vq,th)] if Vq else []
    recs,zp,_=run_beam(nm,src,N=N,B=Bf,I=I*1e-3,zplanes=zp,seed=1,beam_kw=dict(dR=0.0),Efield=EL,iters=int(os.environ.get("ITERS","5")))
    b=best(recs,zp,amax=0.025)
    w=abs(zp-2.415)<=0.05; s=stage({k:v[w] for k,v in recs.items()},zp[w],0.021)
    return dict(I=I,b2=b2,a2=a2,Vq=Vq,th=th,T=b['T'],rms=b['rms87'],z=b['z'],sT=s['T'])
if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='scan0':   # zero current, b2=1e-3 at 104 mm
        J=[(0,'srcB_6mA_V17400.pkl',1.89,1e-3,0,0,0,800)]+[(0,'srcB_6mA_V17400.pkl',1.89,1e-3,0,V,np.deg2rad(t),800) for t in (0,22.5,45,67.5,90,112.5,135,157.5) for V in (100,)]
    else:
        J=[tuple(json.loads(x)) for x in sys.argv[2:]]
    with Pool(4) as p:
        for r in p.imap(job,J): print('I=%g b2=%g a2=%g Vq=%g th=%.1f  rms87=%.2f mm best T=%.3f z=%.3f stageT(21mm)=%.3f'%(r['I'],r['b2'],r['a2'],r['Vq'],np.rad2deg(r['th']),r['rms']*1e3,r['T'],r['z'],r['sT']),flush=True)
