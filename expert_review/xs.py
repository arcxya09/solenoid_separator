import numpy as np
from scipy.integrate import quad
a0=0.529177e-10; k=14.3996e-10  # e^2/4pi eps0 in eV*m
def zbl(Z1,Z2):
    a=0.8854*a0/(Z1**0.23+Z2**0.23); c=[0.1818,0.5099,0.2802,0.02817]; d=[3.2,0.9423,0.4029,0.2016]
    def dV(r):  # dV/dr in eV/m
        x=r/a; phi=sum(ci*np.exp(-di*x) for ci,di in zip(c,d)); dphi=sum(-ci*di/a*np.exp(-di*x) for ci,di in zip(c,d))
        return Z1*Z2*k*(-phi/r**2+dphi/r)
    return dV
def chi(b,Ec,dV):  # impulse approximation CM deflection
    f=lambda u: dV(np.sqrt(b*b+u*u))*b/np.sqrt(b*b+u*u)
    I=2*quad(f,0,50*b+2e-9,limit=200)[0]
    return -I*b/(2*Ec)/b*1.0  # chi = -(1/(2Ec)) * int dV/dr * (b/r) dz  (both halves)
def sig_gt(theta_lab,E,m1,m2,Z1,Z2):
    dV=zbl(Z1,Z2); Ec=E*m2/(m1+m2)
    bs=np.logspace(-11,-9,300); ch=np.minimum(np.array([chi(b,Ec,dV) for b in bs]),1.0)
    th=np.arctan2(np.sin(ch),np.cos(ch)+m1/m2)
    # th decreasing with b
    return np.array([np.pi*np.interp(t,th[::-1],bs[::-1])**2 for t in theta_lab])
TH=np.array([2,5,10,20,50,100])*1e-3
P85=dict(zip([2,5,10,20,50,100],[0,0,0.0085,0.0652,0.0693,0.0693*0.25]))
if __name__=='__main__':
    kB=1.380649e-23
    for gas,m2,Z2,mult in (('He',4,2,1),('N2',14,7,2),('Ar',40,18,1),('Xe',131,54,1)):
        for E in (10e3,30e3):
            s=sig_gt(TH,E,86.9,m2,37,Z2)*mult
            # contamination per Pa: integrate p_hole over angle bins (piecewise const average of bin edges)
            pb=[0.5*(P85[a]+P85[b]) for a,b in zip([2,5,10,20,50],[5,10,20,50,100])]
            sb=s[:-1]-s[1:]; nL=1/(kB*300)*1.97
            leak=nL*np.sum(sb*pb)+nL*s[-1]*P85[100]
            print('%-3s %2d keV sigma(>θ) m^2:'%(gas,E/1e3),' '.join('%.1e'%x for x in s),
                  '| P(>10mrad) at 2.7e-3 Pa = %.2e | Rb85 into hole per ion per Pa = %.2e -> abundance loss at 2.7e-3 Pa: %.2e'%(s[2]*nL*2.7e-3,leak,2.59*leak*2.7e-3))
