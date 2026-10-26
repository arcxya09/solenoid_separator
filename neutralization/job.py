import sys,json,pickle,time
from nz import *
import nz
N=6000
def tube_avg(n):
    """spread n uniformly along each flux tube over the z-range where the tube is threaded (slow ions)"""
    nb=len(BINS)-1; tot=np.zeros(nb); cnt=np.zeros(nb)
    thr=(n>1e-3*n.max())
    w=(2*np.pi*np.maximum(rg,DR/4))[:,None]*np.ones_like(n)  # volume weights
    np.add.at(tot,TB.ravel(),(n*w).ravel()); 
    # volume of tube occupied: cells of the tube between first and last threaded z
    occ=np.zeros_like(n,bool)
    for b in np.unique(TB[thr]):
        m=TB==b; zz=np.where((m&thr).any(0))[0]; occ[:,zz.min():zz.max()+1]|=m[:,zz.min():zz.max()+1]
    np.add.at(cnt,TB.ravel(),(occ*w).ravel())
    return np.where(occ,(tot/np.maximum(cnt,1e-30))[TB],0)
def job(cfg):
    t=time.time(); I=cfg['I']; Te=cfg.get('Te',3.0); Ph=cfg.get('Phi',10.0); eta=cfg.get('eta',0.0); mode=cfg.get('mode','boltz')
    bm=beams(N); recs=run_tracks(bm,None)
    ni=deposit(recs,[I*AB['85']/N,I*AB['87']/N]); phi=None; hist=[]
    for it in range(cfg.get('iters',4)):
        ntot=ni+eta*tube_avg(ni) if eta>0 else ni
        if mode=='boltz':
            phi,ne,o=solve(ntot,Te,Ph,phi0=phi)
        else:  # unneutralized, or uniform fraction f
            phi=sla.spsolve(L.tocsc(),-(qe/eps0)*(1-cfg.get('f',0))*ntot.ravel()*INNER).reshape(NR,NZ)
        recs=run_tracks(bm,efield(phi)); m=evaluate(recs); hist.append(m['z'])
        nnew=deposit(recs,[I*AB['85']/N,I*AB['87']/N]); ni=0.5*ni+0.5*nnew
    jc=np.argmin(abs(zg-1.0)); b=ni[:,jc]>0.05*ni[:,jc].max()
    m.update(cfg); m.update(zhist=hist,err=nz.ERR if mode=='boltz' else 0,phi_coil_beam=[float(phi[b,jc].min()),float(phi[b,jc].max())],phimax=float(phi.max()),phimin=float(phi.min()),t=time.time()-t)
    print(json.dumps(m),flush=True); return m
if __name__=='__main__':
    from multiprocessing import Pool
    cfgs=json.loads(open(sys.argv[1]).read())
    with Pool(int(sys.argv[3]) if len(sys.argv)>3 else 4) as p: res=p.map(job,cfgs,chunksize=1)
    json.dump(res,open(sys.argv[2],'w'),indent=1)
