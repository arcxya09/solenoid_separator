import numpy as np
from sep import loop_field
k=np.sqrt(3); L=1.0*k; a_in=0.13*k; a_out=0.16*k; nz=433; nr=4; B0=1.912
zc=(np.arange(nz)+0.5)/nz*L; ac=a_in+(np.arange(nr)+0.5)/nr*(a_out-a_in)
def sol(r,z):
    br=0;bz=0
    for z0 in zc:
        for a in ac: x,y=loop_field(a,z0,r,z); br=br+x; bz=bz+y
    return br,bz
_,c=sol(np.array([0.0]),np.array([L/2])); S=B0/c[0]
z=np.linspace(-1.31,-1.05,53); rr=np.array([0.0,0.06,0.09,0.104,0.12])
R,Z=np.meshgrid(rr,z,indexing='ij'); Br,Bz=sol(R.ravel(),Z.ravel()); Br=Br.reshape(R.shape)*S; Bz=Bz.reshape(R.shape)*S
print('solenoid only: Bz at r=104mm from %.1f to %.1f mT; |Br| max %.1f mT'%(Bz[3].min()*1e3,Bz[3].max()*1e3,abs(Br[1:]).max()*1e3))
def coils(params):  # list of (a, z0, NI)
    br=np.zeros_like(R); bz=np.zeros_like(R)
    for a,z0,NI in params:
        x,y=loop_field(a,z0,R.ravel(),Z.ravel()); br+=x.reshape(R.shape)*NI; bz+=y.reshape(R.shape)*NI
    return br,bz
def fit(geom):  # geom: list of (a,z0); least squares NI to null Bz at r<=0.12
    cols=[]
    for a,z0 in geom:
        _,y=coils([(a,z0,1.0)]); cols.append(y.ravel())
    A=np.array(cols).T; NI=np.linalg.lstsq(A,-Bz.ravel(),rcond=None)[0]
    br,bz=coils([(a,z0,n) for (a,z0),n in zip(geom,NI)])
    tz=Bz+bz; tr=Br+br
    return NI,abs(tz).max(),abs(tr).max(),tz,tr
for name,geom in (('single loop a=0.20 at z=-1.18',[(0.20,-1.18)]),
                  ('single loop a=0.20 at emitter z=-1.294',[(0.20,-1.294)]),
                  ('pair a=0.20 at -1.33/-1.03',[(0.20,-1.33),(0.20,-1.03)]),
                  ('pair a=0.25 at -1.36/-1.00',[(0.25,-1.36),(0.25,-1.00)]),
                  ('triple a=0.25 at -1.40/-1.18/-0.96',[(0.25,-1.40),(0.25,-1.18),(0.25,-0.96)])):
    NI,mz,mr,tz,tr=fit(geom)
    print('%-38s NI='%name+','.join('%.0f'%n for n in NI)+' A-turn | max|Bz| %.2f mT, max|Br| %.2f mT (r<=120mm, z -1.31..-1.05)'%(mz*1e3,mr*1e3))
    if name.startswith('single loop a=0.20 at emitter'):
        print('   Bz at r=104 along z (mT):',np.round(tz[3][::8]*1e3,1))
