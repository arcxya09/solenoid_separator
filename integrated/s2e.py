# start-to-end job runner: python3 s2e.py jobs.json out.json [nproc]
# job keys: design, src, N, B0, I (total A), f0 (neutralisation), err (ErrField kw), beam_kw (source errors),
#           steer [(thx,thy),(thx,thy)], collector (bool), Isep (pocket SC current), fpocket, env (bool), tag, seed
from collector import *
import sys
def job(c):
    c=dict(c); d=DESIGNS[c['design']]; design=c['design']
    sol=get_sol(design,c.get('B0'))
    B=ErrField(sol,**c['err']) if c.get('err') else sol
    Ef=[steerers(design,c['steer'])] if c.get('steer') else []
    f0=c.get('f0',0.0); I=c.get('I',0.0)
    fneu=fneu_profile(design,f0,z_on=c.get('f_zon'),z_off=c.get('f_zoff'),ramp=c.get('f_ramp',0.02)) if f0 else None
    zp=np.arange(d['zc']-0.30,d['zc']+0.20,0.0025)
    recs,zp,scd=run_beam(design,c['src'],N=c.get('N',4000),B=B,B0=c.get('B0'),I=I,fneu=fneu,Efield=Ef,
                         seed=c.get('seed',0),zplanes=zp,beam_kw=c.get('beam_kw'),iters=c.get('iters',5))
    res=dict(tag=c.get('tag',''),cfg={k:v for k,v in c.items() if k not in ('src',)},src=os.path.basename(c['src']))
    res['best']=best(recs,zp); res['fixed']=fixed(recs,zp,d['zc'],d['a'])
    w=abs(zp-d['zc'])<=0.05                     # collector z stage +-50 mm, x/y follows the Rb-87 centroid
    res['stage']=stage({k:v[w] for k,v in recs.items()},zp[w],d['a'])
    # focus position and spot size vs z (Rb-87), Rb-85 inner radius
    r7=np.hypot(recs['87'][...,0],recs['87'][...,1]); r5=np.hypot(recs['85'][...,0],recs['85'][...,1])
    res['zscan']=dict(z=zp[::2].round(4).tolist(),T=np.nanmean(r7<d['a'],1)[::2].round(4).tolist(),
                      L=np.mean(r5<d['a'],1)[::2].round(5).tolist(),rms87=np.sqrt(np.nanmean(r7**2,1))[::2].round(5).tolist())
    res['sc']=dict(conv=scd.get('conv'),phimax=scd.get('phimin'))
    res['loss']={nm:float(np.mean(~np.isfinite(recs[nm][-1][:,0]))) for nm in recs}
    if c.get('collector'):
        Efc=list(Ef)+([scd['Ef']] if scd.get('Ef') else [])
        col,f7,f5=run_collector(design,recs,zp,B,Isep=c.get('Isep',0.0),fpocket=c.get('fpocket',0.0),Efield=Efc,gk=c.get('colgeom'))
        res['collector']=summarize(col,f7,f5,recs['87'].shape[1],recs['85'].shape[1]); res['collector']['P']={k:(v if isinstance(v,str) else float(v)) for k,v in col.P.items()}
    if c.get('env'):
        # full-line envelope: rerun with dense planes and no SC change (uses the converged SC field)
        ze=np.arange(d['zs'],d['zc']+0.1,0.01)
        sol2=sol; EL=list(Ef)+([scd['Ef']] if scd.get('Ef') else [])
        env={}
        for nm in ('85','87'):
            X,V,qm=beam(sol2,design,nm,1500,c['src'],seed=5,**(c.get('beam_kw') or {}))
            rec=track(B,X,V,qm,ze,ion_dt(design,qm,sol2.B0),Efield=EL,rpipe=d['rpipe'])
            rr=np.hypot(rec[...,0],rec[...,1])
            env[nm]=np.nanpercentile(rr,[5,50,95],axis=1).T.round(5).tolist()
        res['env']=dict(z=ze.round(3).tolist(),**env,Bz=np.interp(ze,sol.zg,sol.Bz[0]).round(5).tolist())
    return res
def show(r):
    b=r['best']; f=r['fixed']; s=r['stage']
    txt='%-28s best T=%.4f a=%.2f z=%.4f rms=%.2f r87max=%.2f r85min=%.2f | fixed T=%.4f pur=%.5f | stage T=%.4f z=%.3f'%(
        r['tag'],b['T'],b['a']*1e3,b['z'],b['rms87']*1e3,b['r87max']*1e3,b['r85min']*1e3,f['T'],f['pur'],s['T'],s['z'])
    if 'collector' in r: txt+=' | col87 '+str({k:(round(v,4) if isinstance(v,float) else v) for k,v in r['collector']['87'].items()})
    return txt
if __name__=='__main__':
    from multiprocessing import Pool
    cf=json.load(open(sys.argv[1])); np_=int(sys.argv[3]) if len(sys.argv)>3 else 4
    with Pool(np_) as p: res=p.map(job,cf,chunksize=1)
    json.dump(res,open(sys.argv[2],'w'))
    for r in res: print(show(r))
