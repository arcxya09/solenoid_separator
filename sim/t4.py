from parax import *
E=10e3; L=1.0
for B0 in np.arange(2.0,4.01,0.05):
    m87,a87=tmap(B0,M87,E,-0.6,L+0.3); m85,a85=tmap(B0,M85,E,-0.6,L+0.3)
    Ds=np.linspace(1.05,2.5,800)
    f=np.array([ (m87[0,0]+(d-1.3)*m87[1,0]) for d in Ds])
    sg=np.where(np.diff(np.sign(f)))[0]
    for i in sg[:1]:
        d=Ds[i]; Dr=np.array([[1,d-1.3],[0,1]]); M7=Dr@m87; M5=Dr@m85
        print(f"B={B0:.2f} a87={a87:.3f} zdet={d:.3f} M21_87={M7[1,0]:.2f} M11_85={M5[0,0]:.4f} sep(R=30mm)={0.03*M5[0,0]*1e3:.2f}mm M12={M7[0,1]:.3f}")
