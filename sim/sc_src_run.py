import os,subprocess,json,sys
from multiprocessing import Pool
def j(a):
    B,I=a; o=subprocess.run(['python3','sc_src.py',str(I)],env=dict(os.environ,B0=str(B)),capture_output=True,text=True)
    try: return json.loads(o.stdout.strip().splitlines()[-1])
    except Exception: return dict(err=o.stderr[-500:],B0=B,I=I)
cf=json.loads(sys.argv[1])
with Pool(4) as p: res=p.map(j,cf)
old=json.load(open('sc_src.json')) if os.path.exists('sc_src.json') else []
json.dump(old+res,open('sc_src.json','w'))
for r in res: print({k:r[k] for k in r if k in ('B0','I','T_0.999','z','rms87','err')})
