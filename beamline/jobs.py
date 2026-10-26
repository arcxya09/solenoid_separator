# generic job runner: python3 jobs.py cfgs.json out.json [nproc]
from bl import *
import json, sys
from multiprocessing import Pool
COL=json.load(open('collector.json')) if os.path.exists('collector.json') else None
def job(cfg):
    cfg=dict(cfg); tag=cfg.pop('tag',''); B0=cfg.pop('B0',2.78); N=cfg.pop('N',6000)
    mis=cfg.pop('mis',None); st=cfg.pop('steer',None)
    sol=get_sol(B0); sl=Mis(sol,**mis) if mis else sol
    out=run(sl,N=N,Ef=steer_field(st) if st else None,**cfg)
    res=dict(tag=tag,B0=B0,mis=mis,steer=st,**cfg); res['best']=best(out)
    if COL:
        res['fixed']=fixed(out,COL['z'],COL['a'])
        res['fixed_rc']=fixed(out,COL['z'],COL['a'],cx=res['best']['cx'],cy=res['best']['cy'])
    return res
if __name__=='__main__':
    cf=json.load(open(sys.argv[1])); np_=int(sys.argv[3]) if len(sys.argv)>3 else 4
    with Pool(np_) as p: res=p.map(job,cf)
    json.dump(res,open(sys.argv[2],'w'),indent=1)
    for r in res:
        b=r['best']; f=r.get('fixed',{})
        print(r['tag'],'| best T=%.4f a=%.2fmm z=%.3f rms=%.2fmm c=(%.2f,%.2f)mm'%(b['T'],b['a']*1e3,b['z'],b['rms87']*1e3,b['cx']*1e3,b['cy']*1e3),
              '| fixed T=%.4f pur=%.5f'%(f.get('T',0),f.get('pur',0)) if f else '')
