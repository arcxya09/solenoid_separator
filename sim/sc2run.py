import os,sys,json,subprocess
from multiprocessing import Pool
def j(a):
    B,I=a
    out=subprocess.run(['python3','-c',f'import sc2; import json; print(json.dumps(sc2.run({I})))'],env=dict(os.environ,B0=str(B)),capture_output=True,text=True)
    d=json.loads(out.stdout.strip().splitlines()[-1]); d['B0']=B; return d
cf=[(B,I) for I in (3e-4,1e-3) for B in (2.80,2.86,2.92,2.98)]
with Pool(4) as p: res=p.map(j,cf)
json.dump(res,open('sc2.json','w'),indent=1)
for r in res: print(r)
