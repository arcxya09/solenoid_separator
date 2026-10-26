from core import *
s=get_sol(1.0); zg=s.zg; b1=s.Bz0[0]   # on-axis field per 1T central
def p_of(M,E): return np.sqrt(2*M*amu*E*qe)
def tmap(B0,M,E,z1,z2):
    """Larmor-frame paraxial 2x2 map from z1 to z2 and Larmor rotation angle"""
    K=qe*B0*b1/(2*p_of(M,E)); dz=zg[1]-zg[0]
    m=np.eye(2); alpha=0
    zz=np.arange(z1,z2,dz)
    Ks=np.interp(zz+dz/2,zg,K,left=0,right=0)
    for k in Ks:
        if k<1e-9: e=np.array([[1,dz],[0,1]])
        else: c,sn=np.cos(k*dz),np.sin(k*dz); e=np.array([[c,sn/k],[-k*sn,c]])
        m=e@m; alpha+=k*dz
    return m,alpha
