import numpy as np, json
kB=1.380649e-23; qe=1.602e-19
def PRb(T): return 10**(9.863-4215/T) if T<312.45 else 10**(9.318-4040/T)   # Pa, solid/liquid (CRC fit)
L=1.97   # source exit -> collector path, m
sig_cx_Rb=3e-18  # m^2 resonant Rb+ + Rb at 10 keV (order 3e-14 cm^2)
sig_cx_gas=1e-19 # m^2 Rb+ + N2/H2O charge exchange (order 1e-15 cm^2)
sig_sc5=3e-20    # m^2 elastic scattering >5 mrad on N2 (screened Coulomb estimate)
out={'T':[], 'P':[], 'nsl':[]}
for T in [200,220,230,240,250,260,273,298]:
    P=PRb(T); n=P/(kB*T); out['T'].append(T); out['P'].append(P); out['nsl'].append(n*sig_cx_Rb*L)
    print(T,'K  P_Rb=%.2e Pa  n sigma L=%.2e'%(P,n*sig_cx_Rb*L))
for P in [1e-4,1e-5,1e-6]:
    n=P/(kB*300); print('gas',P,'Pa  cx %.1e  scat>5mrad %.1e'%(n*sig_cx_gas*L,n*sig_sc5*L))
# beam / source fluxes
I=0.3e-3; Nion=I/qe; print('ions/s %.2e, Rb mass g/day %.3f'%(Nion,Nion*85.5*1.66e-24*86400))
json.dump(out,open('vac.json','w'))
