import numpy as np
qe=1.602e-19; me=9.109e-31; eps0=8.854e-12; kB=1.381e-23; amu=1.6605e-27
vb=np.sqrt(2*30e3*qe/(87*amu)); print('v_beam %.3g m/s, equiv electron energy %.2f eV'%(vb,0.5*me*vb**2/qe))
# beam density in coil (ring R~30-50mm converging, width ~2mm) and drift region
for I in [3e-3,10e-3,30e-3,100e-3]:
    lam=I/vb; n=lam/qe/(2*np.pi*0.04*2e-3)
    phi0=lam/(2*np.pi*eps0)*np.log(0.13/0.04)
    for Te in [1,3,5]:
        lD=np.sqrt(eps0*Te/(n*qe))
    wD=n*qe/(2*eps0*2.8); wD_drift=n*qe/(2*eps0*0.05)
    print('I=%3d mA n=%.2e m^-3 phi0(unneutral, Gauss)=%.0f V lD(3eV)=%.2f mm  w_diocotron 2.8T %.1e /s, 50mT %.1e /s'%(I*1e3,n,phi0,np.sqrt(eps0*3/(n*qe))*1e3,wD,wD_drift))
# gas: rates per beam ion
for P in [1e-5,1e-4,1e-3,2.7e-3]:
    ng=P/(kB*300)
    for se,sc in [(1e-21,1e-21),(1e-20,1e-20)]:
        nu_e=ng*se*vb; nu_s=ng*(se+sc)*vb
        print('P=%.0e Pa ng=%.1e sig_e=%.0e: nu_e=%.2e /s (fill %.2g ms), slow-ion prod/beam ion over 2.5m %.1e, CX loss per pass(sig 1e-19 conservative) %.1e'%(P,ng,se,nu_e,1e3/nu_e,ng*(se+sc)*2.5,ng*1e-19*2.5))
# electron gyroradius and collisions
for B in [0.02,0.05,0.3,2.8]:
    print('B=%.2f T: r_L(e,5eV)=%.3g mm  r_L(N2+,5eV)=%.3g mm r_L(Rb+,5eV)=%.3g mm'%(B,3.37e-6*np.sqrt(5)/B*1e3,np.sqrt(2*28*amu*5*qe)/(qe*B)*1e3,np.sqrt(2*85*amu*5*qe)/(qe*B)*1e3))
n=1e15; Te=3; lnL=10
nu_ee=2.91e-12*n*1e-6*lnL*Te**-1.5*1e6/1e6  # NRL: 2.91e-6 n[cm^-3] lnL T^-1.5
nu_ee=2.91e-6*(n*1e-6)*lnL*Te**-1.5
print('nu_ee(n=1e15,3eV)=%.2e /s ; bounce rate (2.5 m, 3 eV)=%.2e /s'%(nu_ee,np.sqrt(2*Te*qe/me)/5))
# Pastukhov-like: needed confinement tau_e=1/nu_e; solve x e^x ~ nu_ee/nu_e*... Delta/Te
for nu_e in [10,100,1000]:
    from scipy.optimize import brentq
    x=brentq(lambda x: np.sqrt(np.pi)/2*x*np.exp(x)-nu_ee/nu_e,0.01,50)
    print('nu_e=%g: barrier Delta ~ %.1f Te'%(nu_e,x))
# slow ion dwell: N2+ / Rb+ falling through ~Te along 1.25 m (half length)
for M in [28,85]:
    v=np.sqrt(2*3*qe/(M*amu)); print('slow ion M=%d dwell ~%.2f ms (3 eV), ratio to beam transit %.0f'%(M,1.25/v*1e3,(1.25/v)/(2.5/vb)))
