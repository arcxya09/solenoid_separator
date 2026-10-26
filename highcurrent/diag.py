import os,subprocess,json
from multiprocessing import Pool
code='''
import sc3,fromsrc,numpy as np,json
sc3.ring_beam=lambda s,N,M,E,**kw: fromsrc.src_beam(s,M,"src_v2.pkl",N,rng=np.random.default_rng(int(M)))
print(json.dumps(sc3.run(%g)))'''
def j(a):
    lab,I,zlo,zhi=a
    o=subprocess.run(['python3','-c',code%I],env=dict(os.environ,B0='2.78',ZLO=str(zlo),ZHI=str(zhi)),capture_output=True,text=True)
    try: d=json.loads(o.stdout.strip().splitlines()[-1])
    except Exception: return dict(lab=lab,err=o.stderr[-400:])
    d['lab']=lab; return d
cf=[('全程 1mA',1e-3,-9,9),('仅源到线圈入口 z<0',1e-3,-9,0.0),('仅线圈内 0<z<1.05',1e-3,0.0,1.05),('仅出口后 z>1.05',1e-3,1.05,9)]
with Pool(4) as p: res=p.map(j,cf)
json.dump(res,open('diag.json','w'))
for r in res: print(r.get('lab'),r.get('T_0.999'),r.get('z'),r.get('rms87'),r.get('err'))
