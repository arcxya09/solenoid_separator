# machine B 6 mA: Rb-85 impact map on cone/dump, power density, and sputtered-neutral MC into the product hole
from collector import *
import neut2
d=DESIGNS['mB']; sol=get_sol('mB')
recs,zp,scd=run_beam('mB','src/mB_6mA_Vl17.4.pkl',N=6000,I=6e-3,iters=6,zplanes=np.arange(d['zc']-0.3,d['zc']+0.02,0.0025))
col,f7,f5=run_collector('mB',recs,zp,sol,Isep=6e-3,gk=dict(style='pocket',Vland=3000.,h=0.5e-3,tout_deg=float(os.environ.get('TOUT','35'))),Efield=[scd['Ef']],iters=4)
o=summarize(col,f7,f5,6000,6000); print(json.dumps(o['85'])); print(json.dumps(o['87']))
np.savez(os.environ.get('F5OUT','mB_f5.npz'),kc=col.P['tout'],X=f5['X'],V=f5['V'],hit=f5['hit'],names=np.array(col.names),zc=col.zc,a=col.P['a'],Lc=col.P['Lc'],Rd=col.P['Rd'],zE2=col.P['zE2'])
if os.environ.get('SAVEONLY'): raise SystemExit
nm=col.names; X=f5['X']; r=np.hypot(X[:,0],X[:,1]); zz=X[:,2]-col.zc; lab=f5['hit']
I85=6e-3*AB['85']; P=I85*30e3; w=P/6000          # W per macro-ion (all Rb-85 at ~30 keV)
out=dict(frac={nm[l] if l>0 else 'upstream':float(np.mean(lab==l)) for l in set(lab.tolist())})
for part in ('cone','dump','housing','E2'):
    m=lab==nm.index(part)
    if m.any(): out[part+'_r_mm']=np.percentile(r[m]*1e3,[0,5,50,95,100]).round(1).tolist()
md=lab==nm.index('dump'); edges=np.linspace(0,0.2,41); h,_=np.histogram(r[md],edges)
area=np.pi*(edges[1:]**2-edges[:-1]**2); pd=h*w/area/1e4    # W/cm^2 (azimuthally averaged)
out['dump_W']=float(md.sum()*w); out['cone_W']=float(np.sum(lab==nm.index('cone'))*w)
out['dump_Wcm2_profile']=dict(r_mm=(edges[:-1]*1e3).round(0).tolist(),W_cm2=pd.round(3).tolist())
mc=lab==nm.index('cone')
if mc.any():
    # cone outer surface area element: 2 pi r dl, dl=dz*sqrt(1+kc^2)
    kc=col.P['tout']; zb=np.linspace(0,col.P['Lc'],11); hc,_=np.histogram(zz[mc],zb)
    rr=col.P['a']+0.4e-3+kc*0.5*(zb[1:]+zb[:-1]); ac=2*np.pi*rr*np.diff(zb)*np.sqrt(1+kc**2)
    out['cone_Wcm2_peak']=float((hc*w/ac/1e4).max())
print(json.dumps(out))
# sputtered neutrals: start points = Rb-85 impacts on cone and dump
sel=mc|md
geom=dict(a=col.P['a'],r0=col.P['a']+0.4e-3,kc=float(col.P['tout']),Lc=col.P['Lc'],Rd=col.P['Rd'],Rl=d['rpipe'],
          zlo=-1.0,zE2=col.P['zE2'],rE2=d['rE2'],rE2o=d['rE2']+10e-3)
neu={}
for s in (0.5,0.9,0.99,0.999):
    f,f1=neut2.run(r[sel],np.where(mc[sel],0.0,col.P['Lc']),geom,s,N=400000); neu[str(s)]=dict(total=f,direct=f1); print('stick',s,f,f1)
out['neutral']=neu; out['geom']={k:float(v) for k,v in geom.items()}
json.dump(dict(out=out,col85=o['85'],col87=o['87']),open('mB_dump.json','w'),indent=1)
