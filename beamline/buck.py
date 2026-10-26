# fringe flux at the emitter and a bucking coil that nulls it
from core import *
import json
s=get_sol(2.78)
def flux(sol_br_bz_fun,z,R=0.03,n=60):
    r=(np.arange(n)+0.5)/n*R; _,bz=sol_br_bz_fun(r,np.full(n,z)); return float(np.sum(bz*2*np.pi*r*R/n))
def solf(r,z):
    br,by,bz=s.field(r,0*r,z); return br,bz
res={'zs':[], 'B':[], 'Beff':[]}
for ze in np.arange(-1.0,-0.25,0.05):
    ph=flux(solf,ze); res['zs'].append(round(float(ze),3)); res['Beff'].append(ph/(np.pi*0.03**2)); res['B'].append(float(np.interp(ze,s.zg,s.Bz[0])))
# bucking coil: single coil (radius ab, length 40 mm, 4 layers) centred on the emitter plane, z_e=-0.685
ze=-0.685
def coil(ab,zc,r,z,len_=0.04,nz=8):
    br=0;bz=0
    for zz in zc+(np.arange(nz)+0.5)/nz*len_-len_/2:
        a,b=loop_field(ab,zz,r,z); br=br+a/nz; bz=bz+b/nz
    return br,bz
ph0=flux(solf,ze); out={}
for ab in (0.08,0.12,0.16):
    ph1=flux(lambda r,z: coil(ab,ze,r,z),ze)   # per ampere-turn
    NI=-ph0/ph1
    # residual flux error across the emitter depth (+-2 mm) and field it puts on the solenoid centre
    dz=[flux(lambda r,z: (solf(r,z)[0]+NI*coil(ab,ze,r,z)[0], solf(r,z)[1]+NI*coil(ab,ze,r,z)[1]),ze+d)/(np.pi*0.03**2) for d in (-2e-3,2e-3)]
    _,bc=coil(ab,ze,np.array([0.]),np.array([0.5])); 
    zz=np.linspace(-0.9,-0.4,51); _,bb=coil(ab,ze,np.zeros_like(zz),zz)
    out[str(ab)]=dict(NI=float(NI),Beff_pm2mm=[float(x) for x in dz],B_at_centre=float(NI*bc[0]),z=zz.round(3).tolist(),Btot=(np.interp(zz,s.zg,s.Bz[0])+NI*bb).round(5).tolist())
json.dump(dict(scan=res,ze=ze,phi0=ph0,Beff0=ph0/(np.pi*0.03**2),buck=out),open('buck.json','w'),indent=0)
print('Beff at emitter',ph0/(np.pi*0.03**2))
for k,v in out.items(): print(k,v['NI'],v['Beff_pm2mm'],v['B_at_centre'])
for z,b,be in zip(res['zs'],res['B'],res['Beff']): print(z,round(b*1e3,2),round(be*1e3,2))
