import numpy as np, json
kB=8.617333e-5; qe=1.602176634e-19; amu=1.66053907e-27; kBJ=1.380649e-23
IP=4.177; g=0.5
T=np.arange(900,1801,20.0)
WF={'W（多晶 4.55 eV）':4.55,'Re（4.96 eV）':4.96,'Ir（5.27 eV）':5.27}
ion={}
for k,phi in WF.items():
    a=g*np.exp((phi-IP)/(kB*T)); ion[k]=(a/(1+a)).tolist()
# coverage: ion desorption rate k=nu*exp(-Q/kT); n_ads=J/k ; monolayer 5e14 /cm2
nu=1e13; ML=5e14
J_ours=0.16e-3/qe     # 0.3 mA over 1.9 cm2 -> ions/cm2/s (0.16 mA/cm2)
J_thr=10e-3/qe        # ion-thruster class 10 mA/cm2
cov={}
for Q in (1.8,2.2,2.6):
    k=nu*np.exp(-Q/(kB*T))
    cov['Q=%.1f eV, 0.16 mA/cm²'%Q]=(J_ours/k/ML).tolist()
    cov['Q=%.1f eV, 10 mA/cm²'%Q]=(J_thr/k/ML).tolist()
# Rb supply through porous W: Knudsen
def Prb(Tr): return 101325*10**(4.312-4040/Tr)   # liquid Rb, Pa (Alcock)
m=85.47*amu
Tr=np.arange(330,560,5.0)
Phi_in=Prb(Tr)/np.sqrt(2*np.pi*m*kBJ*Tr)/1e4    # atoms/cm2/s hitting pores
A=1.885   # cm2 emitting area
sup={}
for lab,f in (('透过率 1e-4',1e-4),('透过率 4e-4（典型）',4e-4),('透过率 1.5e-3',1.5e-3)):
    sup[lab]=(Phi_in*f*A*qe*1e3).tolist()   # equivalent ion current mA (beta~1)
# required reservoir T for 0.3 mA with typical transmission
need={}
for I in (0.1,0.3,1.0,2.5):
    cur=Phi_in*4e-4*A*qe*1e3; need[I]=float(np.interp(I,cur,Tr))
# consumption
def mg_day(I_mA,beta=0.99): return I_mA*1e-3/qe/beta*85.47*amu*86400*1e6
cons={I:mg_day(I) for I in (0.1,0.3,1.0)}
# thermal: radiative balance
sig=5.670374e-8
def Prad(Tk,eps,Aeff): return eps*sig*Tk**4*Aeff
R=0.03
front=2*np.pi*R*4e-3          # 4 mm-wide hot front face
back_sides=2*np.pi*R*(4e-3+2*4e-3)
heat={}
for Tk in (1200,1300,1400,1500):
    p_front=Prad(Tk,0.2,front); p_shield=Prad(Tk,0.2,back_sides)/3   # 2-layer shields ~1/3
    heat[Tk]=dict(front=p_front,shielded=p_shield,cond=8.0,total=p_front+p_shield+8.0)
# neutral efflux & CX in extraction region
cx={}
for k,phi in WF.items():
    a=g*np.exp((phi-IP)/(kB*1400)); b=a/(1+a)
    Jn=J_ours*(1-b)/b                         # neutrals /cm2/s leaving
    vbar=np.sqrt(8*kBJ*1400/(np.pi*m))*100    # cm/s
    n0=4*Jn/vbar                               # cm-3 at surface (half-space)
    sig_cx=1.5e-14                             # cm2, resonant Rb+-Rb at keV (order of magnitude)
    Leff=0.5                                   # cm, density falls ~ over ring width scale
    cx[k]=dict(beta=b,neutral_flux=Jn,n0=n0,P_cx=n0*sig_cx*Leff,neutral_mg_day=Jn*A*m*86400*1e6)
out=dict(T=T.tolist(),ion=ion,cov=cov,Tr=Tr.tolist(),sup=sup,need=need,cons=cons,heat=heat,cx=cx)
json.dump(out,open('feas.json','w'),indent=1)
print('need',need); print('cons',cons); print('heat',heat); print(json.dumps(cx,indent=1))
for k in ion: print(k,[round(ion[k][list(T).index(t)],4) for t in (1200,1400,1600)])
for k in cov: 
    c=np.array(cov[k]); print(k,'T for coverage<1%:',T[np.argmax(c<0.01)] if (c<0.01).any() else None)
