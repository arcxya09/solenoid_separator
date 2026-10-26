import os,sys,json
import sc2, fromsrc, numpy as np
def rb(s,N,M,E,**kw): return fromsrc.src_beam(s,M,'src_0p1mA.pkl',N,rng=np.random.default_rng(int(M)))
sc2.ring_beam=rb
I=float(sys.argv[1]); r=sc2.run(I); r['B0']=sc2.B0; print(json.dumps(r))
