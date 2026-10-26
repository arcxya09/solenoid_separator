from parax import *
E=10e3; S=0.6; L=1.0
print('B at source', np.interp(-S,zg,b1))
for B0 in np.arange(2.0,4.01,0.1):
    m87,a87=tmap(B0,M87,E,-S,L+0.6); m85,a85=tmap(B0,M85,E,-S,L+0.6)
    # find D where M12(87)=0 for D in [0.1,1.5]: M(D)=Drift(D-0.6)*m
    Ds=np.linspace(0.62-0.5,1.6,500); best=None
    f=[ (m87[0,1]+(d-0.6)*m87[1,1]) for d in Ds]
    sg=np.where(np.diff(np.sign(f)))[0]
    for i in sg[:1]:
        d=Ds[i]; M7=np.array([[1,d-0.6],[0,1]])@m87; M5=np.array([[1,d-0.6],[0,1]])@m85
        print(f"B={B0:.1f} a87={a87:.3f} D={d:.3f} M11={M7[0,0]:.3f} M12_85={M5[0,1]:.3f} sep(th=.06)={0.06*M5[0,1]*1e3:.1f}mm  dalpha={a85-a87:.3f}")
