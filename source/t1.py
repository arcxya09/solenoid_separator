import time
from es import *
t=time.time(); g=Geom(); print('build',time.time()-t, g.fix.shape); t=time.time(); g.solve(); print('solve',time.time()-t)
X,V,m=emit(g,2000); t=time.time()
out,alive,_=track(g,X,V,m,zout=0.06); print('track',time.time()-t, alive.mean())
res,_=beam_stats(out,g); print({k:(round(v*1e3,4) if isinstance(v,float) else v) for k,v in res.items()})
