from wp import *
from multiprocessing import Pool
def job(a):
    R0,B=a; nm='r_%g_%g'%(R0,B); mk(nm,R0,B,'C3')
    recs,zp,_=run(nm,N=400,zplanes=np.arange(1.0,1.8,0.0025)); zf,rms=focus(recs,zp); return (R0,B,float(zf))
GUESS={(0.03,1):1.10,(0.03,2):1.95,(0.03,3):2.80,(0.06,1):1.10,(0.06,2):1.88,(0.06,3):2.66,(0.08,1):1.09,(0.08,2):1.84,(0.08,3):2.50}
if __name__=='__main__':
    jobs=[(R0,round(B0+dB,3)) for (R0,n),B0 in GUESS.items() for dB in np.arange(-0.06,0.061,0.02)]
    with Pool(4) as p: res=p.map(job,jobs)
    out={}
    for (R0,n),B0 in GUESS.items():
        pts=sorted([(B,zf) for R,B,zf in res if R==R0 and abs(B-B0)<0.07])
        Bs=np.array([p[0] for p in pts]); zs=np.array([p[1] for p in pts])
        ok=(zs>1.01)&(zs<1.79)
        # zf decreases with B
        Bf=float(np.interp(1.365,zs[ok][::-1],Bs[ok][::-1])); slope=float(np.polyfit(Bs[ok],zs[ok],1)[0])
        out['%g_%d'%(R0,n)]=dict(R0=R0,n=n,B=Bf,dzdB=slope,pts=pts); print(R0,n,round(Bf,4),'dz/dB %.2f m/T'%slope, pts)
    json.dump(out,open('refine.json','w'))
