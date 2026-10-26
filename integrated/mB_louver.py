# machine B 6 mA: planar vs deep-groove vs chevron louver dump; sputtered Rb-85 entering the product hole
# input: mB_f5.npz (Rb-85 impacts from mB_dump.py with SAVEONLY=1). usage: python3 mB_louver.py mB_f5.npz out.json
import sys, json, numpy as np
from multiprocessing import Pool
from neut2 import cosdir
import neut3
f=np.load(sys.argv[1],allow_pickle=True); names=list(f['names']); X=f['X'].copy(); V=f['V']; hit=f['hit']
zc=float(f['zc']); a=float(f['a']); Lc=float(f['Lc']); Rd=float(f['Rd']); zE2=float(f['zE2'])
X[:,2]-=zc; NTOT=len(X); PTOT=6e-3*0.7217*30e3          # W of Rb-85 (all macro-ions, 30 keV)
wion=PTOT/NTOT
kc=float(f['kc']) if 'kc' in f else float(np.tan(np.radians(35))); r0=a+0.4e-3; rc=r0+kc*Lc
base=dict(a=a,r0=r0,kc=kc,Lc=Lc,Rd=Rd,Rl=0.20,zlo=-1.0,zE2=zE2,rE2=0.13,rE2o=0.14)
md=hit==names.index('dump'); mc=hit==names.index('cone')
Xd=X[md]; Vd=V[md]; Xd[:,2]=Lc; rd=np.hypot(Xd[:,0],Xd[:,1])
mion=(Xd[:,0]*Vd[:,0]+Xd[:,1]*Vd[:,1])/rd/Vd[:,2]
th_planar=np.degrees(np.arccos(Vd[:,2]/np.linalg.norm(Vd,axis=1)))
Xc=X[mc]; rcn=np.hypot(Xc[:,0],Xc[:,1])
def yfac(cth,p=1.7,cmin=0.2): return np.maximum(cth,cmin)**(-p)
def ionslope(r):   # local median ion slope dr/dz at the dump plane
    e=np.linspace(rc,Rd,30); c=0.5*(e[1:]+e[:-1]); med=np.array([np.median(mion[(rd>=e[i])&(rd<e[i+1])]) if np.sum((rd>=e[i])&(rd<e[i+1]))>5 else np.nan for i in range(29)])
    ok=np.isfinite(med); return np.interp(r,c[ok],med[ok])
def fins(p,D,chev):
    rt=np.arange(rc+p,Rd-1e-4,p); m=ionslope(rt)
    if chev: m=m+chev*p/D       # extra tilt so that the fin bottom passes the next fin top by chev*p (chev=1: opaque)
    return rt,m
CONF={'planar':dict(D=0),
      'G15':dict(D=15e-3,p=5e-3,chev=0,t=0.5e-3),'G25':dict(D=25e-3,p=5e-3,chev=0,t=0.5e-3),
      'C15':dict(D=15e-3,p=5e-3,chev=1.0,t=0.5e-3),'C25':dict(D=25e-3,p=5e-3,chev=1.0,t=0.5e-3),'G25k':dict(D=25e-3,p=5e-3,chev=0,t=0.1e-3),'G40k':dict(D=40e-3,p=5e-3,chev=0,t=0.1e-3)}
def build(name):
    c=CONF[name]
    if c['D']==0: return neut3.Geo(**base),None
    fs=fins(c['p'],c['D'],c['chev']); return neut3.Geo(**base,D=c['D'],fins=fs,rc=rc),fs
def launch(name,rng,N):
    """sputter start points, directions and weights (relative yield) for N atoms, plus ion-power bookkeeping"""
    g,fs=build(name); c=CONF[name]
    # ions: dump ions (incl. fin tips) + cone ions; atoms are sampled per ion (equal number per macro-ion)
    nd,nc_=len(Xd),len(Xc); info={}
    if c['D']==0:
        Pd=Xd.copy(); Nd=np.tile([0,0,-1.],(nd,1)); cd=np.cos(np.radians(th_planar)); sd=np.full(nd,1)
    else:
        Q,S,Nn,cth=neut3.ion_impacts(g,Xd.copy(),Vd)
        tip=rng.random(nd)<c['t']/c['p']
        Pd=np.where(tip[:,None],Xd,Q); Nd=np.where(tip[:,None],np.array([[0,0,-1.]]),Nn)
        cd=np.where(tip,np.cos(np.radians(th_planar)),cth); sd=np.where(tip,-2,S)
        lab=np.where(sd>=100,100,sd); info['ion_land']={str(k):float(np.mean(lab==k)*len(Xd)/NTOT) for k in np.unique(lab)}
        info['ion_lost']=float(np.mean(S<0))
        fin=(sd>=100); info['fin_theta_deg']=np.percentile(np.degrees(np.arccos(cth[fin])),[5,50,95]).round(1).tolist() if fin.any() else None
        bot=(sd==5); info['bot_theta_deg']=np.percentile(np.degrees(np.arccos(cth[bot])),[5,50,95]).round(1).tolist() if bot.any() else None
        # power density on fins (per fin, 5 z bins) and bottom (5 mm annuli)
        pk=0; finW=[]
        for k in range(len(fs[0])):
            w=sd==100+k
            if not w.any(): finW.append(0.0); continue
            zb=np.linspace(Lc,Lc+c['D'],6); h,_=np.histogram(Pd[w,2],zb); rm=fs[0][k]+fs[1][k]*(0.5*(zb[1:]+zb[:-1])-Lc)
            ar=2*np.pi*rm*np.diff(zb)*np.sqrt(1+fs[1][k]**2); pk=max(pk,(h*wion/ar/1e4).max()); finW.append(float(w.sum()*wion))
        info['fin_peak_Wcm2']=float(pk); info['fin_W_max']=float(max(finW)); info['fins_W']=float(sum(finW))
        rb=np.hypot(Pd[bot,0],Pd[bot,1]); e=np.arange(rc,Rd+5e-3,5e-3); h,_=np.histogram(rb,e)
        info['bot_peak_Wcm2']=float((h*wion/(np.pi*(e[1:]**2-e[:-1]**2))/1e4).max()) if bot.any() else 0.0
        info['bot_W']=float(bot.sum()*wion); info['tip_W']=float(tip.sum()*wion)
        info['n_fins']=len(fs[0]); info['fin_r_top_mm']=[round(fs[0][0]*1e3,1),round(fs[0][-1]*1e3,1)]
        info['fin_tilt_deg']=[round(float(np.degrees(np.arctan(fs[1][0]))),1),round(float(np.degrees(np.arctan(fs[1][-1]))),1)]
        info['Yfac_mean_dump']=float(np.mean(yfac(cd)))
    # cone outer wall ions (same in all cases): incidence from V vs cone normal
    ph=np.arctan2(Xc[:,1],Xc[:,0]); zc_=np.clip((rcn-r0)/kc,0,Lc)
    Pc=np.stack([rcn*np.cos(ph),rcn*np.sin(ph),zc_],1); Ncn=np.stack([np.cos(ph),np.sin(ph),np.full(nc_,-kc)],1); Ncn/=np.linalg.norm(Ncn,axis=1)[:,None]
    Vc=V[mc]/np.linalg.norm(V[mc],axis=1)[:,None]; ccn=np.abs(np.sum(Vc*Ncn,1))
    return g,(Pd,Nd,cd,sd),(Pc,Ncn,ccn),info
def job(args):
    name,s,N,seed=args; rng=np.random.default_rng(seed)
    g,(Pd,Nd,cd,sd),(Pc,Ncn,ccn),info=launch(name,rng,N)
    res={}; tipmask=None
    for part,(P0,N0,c0,lab0),nion in (('dump',(Pd,Nd,cd,sd),len(Pd)),('cone',(Pc,Ncn,ccn,np.full(len(Pc),4)),len(Pc))):
        if nion==0: res[part]=dict(p_escape=0.0,F=0.0,Fy=0.0,n=0,by_surface={}); continue
        n=N*nion//(len(Pd)+len(Pc))+1000; ch=50000; acc={}
        for i0 in range(0,n,ch):
            m=min(ch,n-i0); i=rng.integers(0,nion,m); P=P0[i]+0.0; Dd=cosdir(rng,N0[i],m); P+=1e-7*Dd
            got=neut3.transport(g,P,Dd,None,s,rng=rng); w=yfac(c0[i]); lb=np.where(lab0[i]>=100,100,lab0[i])
            for k in np.unique(lb):
                q=lb==k; a_=acc.setdefault(str(k),[0,0.0,0]); a_[0]+=int(got[q].sum()); a_[1]+=float(w[q&got].sum()); a_[2]+=int(q.sum())
        fi=nion/NTOT; tot=sum(v[0] for v in acc.values())
        res[part]=dict(p_escape=tot/n,F=fi*tot/n,Fy=fi*sum(v[1] for v in acc.values())/n,n=n,
                       by_surface={k:dict(hits=v[0],n=v[2],F=fi*v[0]/n,Fy=fi*v[1]/n) for k,v in acc.items()})
    return name,s,res,info
if __name__=='__main__':
    N=int(__import__('os').environ.get('NN','4000000'))
    sel=__import__('os').environ.get('CONFS'); names_=sel.split(',') if sel else list(CONF)
    jobs=[(nm,s,N,11+i) for i,nm in enumerate(names_) for s in (0.85,0.9,0.95,0.99)]
    out={}
    with Pool(4) as pool:
        for nm,s,res,info in pool.imap_unordered(job,jobs):
            out.setdefault(nm,dict(info=info,conf=CONF[nm]))[str(s)]=res
            F=res['dump']['F']+res['cone']['F']; Fy=res['dump']['Fy']+res['cone']['Fy']
            print(nm,s,'F %.2e (dump %.2e cone %.2e)  yield-weighted %.2e'%(F,res['dump']['F'],res['cone']['F'],Fy),flush=True)
            json.dump(out,open(sys.argv[2],'w'),indent=1)
    for nm in out: print(nm,json.dumps(out[nm]['info']))
