from scrun import *
from multiprocessing import Pool
def j(a):
    I,Vl,fn,M=a
    import es
    return run(I,Vl,N=6000,seed=2,save=fn)
with Pool(2) as p: print(p.map(j,[(1e-4,-15000,'src_0p1mA.pkl',None),(1e-3,-18000,'src_1mA.pkl',None)]))
