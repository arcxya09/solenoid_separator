from scrun import *
EX=dict(rnd=0.0014,gl=0.010,t=0.003,zl=0.032)
r=run(1e-4,-16000,N=6000,seed=2,save='src_v2.pkl',extra=EX); print(r)
