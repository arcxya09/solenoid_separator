from wp import *
from multiprocessing import Pool
REF=json.load(open('refine.json')); I1=json.load(open('i1.json'))
def stage_any(recs,zp,P=0.999,zc=1.365,dz=0.05):
    out=None
    for k in np.where(abs(zp-zc)<=dz+1e-9)[0]:
        b=at_plane(recs,zp,k,P)
        if out is None or b['T']>out['T']: out=b
    return out
def job(a):
    key,I,B=a; d=REF[key]; nm='b%s_%g'%(key,B); mk(nm,d['R0'],B,'C3')
    recs,zp,_=run(nm,N=1500,I=I*1e-3,zplanes=np.arange(1.25,1.5,0.005),seed=1,iters=4)
    s=stage_any(recs,zp); s.update(key=key,I=I,B=B); return s
def predB(key,I):
    d=REF[key]; pts=[(0,1.365)]+[(r['I'],r['z']) for r in I1 if r['key']==key and r['T']>0.9]
    Iz=np.array(pts); sl=np.polyfit(Iz[:,0],Iz[:,1],1)[0] if len(pts)>1 else 0
    dzs=np.interp(I,Iz[:,0],Iz[:,1]) if I<=Iz[-1,0] else Iz[-1,1]+sl*(I-Iz[-1,0])
    return d['B']*(1-(dzs-1.365)/(d['dzdB']*-1)/d['B']*-1) if False else d['B']-(dzs-1.365)/d['dzdB']
if __name__=='__main__':
    jobs=[]
    for key in ('0.03_3','0.06_2','0.08_1'):
        for I in (0.1,0.2,0.3,0.5,1.0,1.5,2.0,3.0):
            Bp=predB(key,I)
            for f in (-0.0075,0,0.0075): jobs.append((key,I,round(Bp*(1+f),4)))
    with Pool(4) as p: res=p.map(job,jobs,chunksize=1)
    json.dump(res,open('bscan.json','w'))
    for r in res: print('%-7s I=%.2f B=%.4f T=%.4f a=%.2f z=%.3f rms=%.2f r85min=%.2f'%(r['key'],r['I'],r['B'],r['T'],r['a']*1e3,r['z'],r['rms87']*1e3,r['r85min']*1e3),flush=True)
