# sc3 + variable drift length, einzel pre-compensation knob, finer radial grid, under-relaxation
from bl import *
import bl, json, sys
from multiprocessing import Pool
rg=np.arange(0,0.09,0.00025)
def build_E(zg,recs,ws):
    lam=np.zeros((len(zg),len(rg)))
    for rec,w in zip(recs,ws):
        r=np.hypot(rec[...,0],rec[...,1]); vz=rec[...,5]
        for i in range(len(zg)):
            ok=np.isfinite(r[i]); idx=np.minimum(np.searchsorted(rg,r[i][ok]),len(rg)-1)
            np.add.at(lam[i],idx,w/vz[i][ok])
    lam=np.cumsum(lam,1); Er=np.zeros_like(lam); Er[:,1:]=lam[:,1:]/(2*np.pi*eps0*rg[1:]); return Er
def mkE(zg,Er):
    def Ef(X):
        r=np.hypot(X[:,0],X[:,1]); z=X[:,2]
        fz=np.clip((z-zg[0])/(zg[1]-zg[0]),0,len(zg)-1.001); fr=np.clip(r/(rg[1]-rg[0]),0,len(rg)-1.001)
        i=fz.astype(int); j=fr.astype(int); u=fz-i; t=fr-j
        e=Er[i,j]*(1-u)*(1-t)+Er[i+1,j]*u*(1-t)+Er[i,j+1]*(1-u)*t+Er[i+1,j+1]*u*t
        e=np.where((z<zg[0])|(z>zg[-1]),0,e); rs=np.where(r>0,r,1)
        return np.stack([e*X[:,0]/rs,e*X[:,1]/rs,0*e],1)
    return Ef
def einzel(X,V,dV):
    """einzel retune by dV kV (more negative): mean r' -0.4 mrad/kV, r'-r slope -2.5 /m per kV (source/vlscan.json)"""
    r=np.hypot(X[:,0],X[:,1]); ur=np.stack([X[:,0]/r,X[:,1]/r],1); dr=r-r.mean()
    drp=(-0.4e-3*dV-2.5*dV*dr)*V[:,2]
    V=V.copy(); V[:,:2]+=drp[:,None]*ur; return V
def run(a):
    I,B0,zs,dV=a[:4]; N=a[4] if len(a)>4 else 3000
    s=get_sol(B0); zg=np.arange(zs,1.9,0.0025)
    beams={}
    for nm in ('85','87'):
        X,V,qm=beam(s,nm,N,src='src_0p3mA.pkl',zs=zs,seed=3); beams[nm]=(X,einzel(X,V,dV) if dV else V,qm)
    Er=None; Ef=None; hist=[]
    NIT=int(__import__('os').environ.get('NIT','6')); RLX=float(__import__('os').environ.get('RLX','0.5')); zh=[]
    for it in range(NIT if I>0 else 1):
        recs={nm:track(s,*beams[nm],zg,0.02/(beams[nm][2][0]*B0),Efield=Ef) for nm in beams}
        if I>0:
            En=build_E(zg,[recs['85'],recs['87']],[I*AB['85']/N,I*AB['87']/N])
            if Er is not None: hist.append(float(np.sqrt(np.mean((En-Er)**2)/np.mean(En**2))))
            sl=zg>=1.0; zh.append(best({nm:recs[nm][sl] for nm in recs},zg[sl])['z'])
            Er=En if Er is None else (1-RLX)*Er+RLX*En; Ef=mkE(zg,Er)
    sel=zg>=1.0; out={nm:recs[nm][sel] for nm in recs}; zp=zg[sel]
    res=dict(I=I,B0=B0,zs=zs,dV=dV,conv=hist,zhist=zh,best=best(out,zp))
    r7=np.hypot(out['87'][...,0],out['87'][...,1]); r5=np.hypot(out['85'][...,0],out['85'][...,1])
    # practical metric: hole radius <=6.5 mm, any z; best T with purity>=0.999 (0/N stats)
    bestp=(0,0,0)
    for k in range(0,len(zp),2):
        for a_ in np.arange(0.004,0.0066,0.0002):
            n7=np.nanmean(r7[k]<a_); n5=np.nanmean(r5[k]<a_)
            if pur(n7,n5)>=0.999 and n7>bestp[0]: bestp=(float(n7),float(zp[k]),float(a_))
    res['practical']=bestp
    res['zscan']=dict(z=zp[::2].round(4).tolist(),T=np.nanmean(r7<0.0059,1)[::2].round(4).tolist(),L=np.nanmean(r5<0.0059,1)[::2].round(5).tolist())
    # ring width at coil entrance
    k0=np.argmin(abs(zg-0.0)); rr=np.hypot(recs['87'][k0][:,0],recs['87'][k0][:,1]); res['w_entr']=float(np.nanpercentile(rr,95)-np.nanpercentile(rr,5))
    return res
if __name__=='__main__':
    cf=json.load(open(sys.argv[1]))
    with Pool(4) as p: res=p.map(run,[tuple(c) for c in cf])
    json.dump(res,open(sys.argv[2],'w'))
    for r in res: print(r['I']*1e3,'mA B',r['B0'],'zs',r['zs'],'dV',r['dV'],'| best T=%.4f z=%.3f a=%.2f'%(r['best']['T'],r['best']['z'],r['best']['a']*1e3),'| practical',r['practical'],'| w95 entr %.1f mm'%(r['w_entr']*1e3),'conv',[round(x,3) for x in r['conv']],'zh',[round(x,3) for x in r['zhist']])
