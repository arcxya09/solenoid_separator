import refine, json, numpy as np
from multiprocessing import Pool
G={(0.045,2):1.93,(0.07,2):1.85}
if __name__=='__main__':
    jobs=[(R0,round(B0+dB,3)) for (R0,n),B0 in G.items() for dB in np.arange(-0.06,0.061,0.02)]
    with Pool(4) as p: res=p.map(refine.job,jobs)
    out=json.load(open('refine.json'))
    for (R0,n),B0 in G.items():
        pts=sorted([(B,zf) for R,B,zf in res if R==R0]); Bs=np.array([p[0] for p in pts]); zs=np.array([p[1] for p in pts]); ok=(zs>1.01)&(zs<1.79)
        Bf=float(np.interp(1.365,zs[ok][::-1],Bs[ok][::-1])); sl=float(np.polyfit(Bs[ok],zs[ok],1)[0])
        out['%g_%d'%(R0,n)]=dict(R0=R0,n=n,B=Bf,dzdB=sl,pts=pts); print(R0,n,Bf,sl,pts)
    json.dump(out,open('refine.json','w'))
