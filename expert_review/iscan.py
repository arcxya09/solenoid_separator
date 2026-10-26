from wp import *
from multiprocessing import Pool
REF=json.load(open('refine.json'))
def metrics(recs,zp):
    b=best(recs,zp)
    k=int(np.argmin(abs(zp-b['z'])))
    r7=np.hypot(recs['87'][k][:,0],recs['87'][k][:,1]); r5=np.hypot(recs['85'][k][:,0],recs['85'][k][:,1])
    r7=r7[np.isfinite(r7)]; r5=r5[np.isfinite(r5)]
    b['r87_99']=float(np.percentile(r7,99)); b['r85_01']=float(np.percentile(r5,1)) if len(r5) else np.nan
    b['loss87']=float(1-np.isfinite(recs['87'][-1][:,0]).mean())
    # fixed collector at design plane with a=5.9 mm hole... use design-plane focus metrics too
    k0=int(np.argmin(abs(zp-1.365)))
    a7=np.hypot(recs['87'][k0][:,0],recs['87'][k0][:,1]); a5=np.hypot(recs['85'][k0][:,0],recs['85'][k0][:,1])
    b['fix_T']=float(np.mean(a7<5.9e-3)); b['fix_L85']=float(np.mean(a5<5.9e-3))
    return b
def job(a):
    key,I=a; d=REF[key]; nm='d'+key; mk(nm,d['R0'],d['B'],'C3')
    recs,zp,_=run(nm,N=3000,I=I*1e-3,zplanes=np.arange(1.0,2.2,0.005),seed=1)
    m=metrics(recs,zp); m.update(key=key,I=I); return m
if __name__=='__main__':
    keys=sys.argv[1].split(','); Is=[float(x) for x in sys.argv[2].split(',')]
    with Pool(4) as p: res=p.map(job,[(k,I) for k in keys for I in Is],chunksize=1)
    json.dump(res,open(sys.argv[3],'w'))
    for r in res: print('%-7s I=%.2f T=%.4f a=%.2f z=%.3f rms=%.2f r87_99=%.2f r85_01=%.2f fixT=%.3f fixL=%.4f'%(r['key'],r['I'],r['T'],r['a']*1e3,r['z'],r['rms87']*1e3,r['r87_99']*1e3,r['r85_01']*1e3,r['fix_T'],r['fix_L85']))
