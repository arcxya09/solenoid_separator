# generic parallel runner: args = list of [script, I, {env}]
import os,subprocess,json,sys
from multiprocessing import Pool
def j(a):
    sc,I,env=a
    o=subprocess.run(['python3',sc,str(I)],env=dict(os.environ,**{k:str(v) for k,v in env.items()}),capture_output=True,text=True)
    try: d=json.loads(o.stdout.strip().splitlines()[-1])
    except Exception: d=dict(err=o.stderr[-400:])
    d.update(script=sc,env=env,Iarg=I); return d
cf=json.loads(sys.argv[1]); out=sys.argv[2]
with Pool(int(os.environ.get('NP','4'))) as p: res=p.map(j,cf)
old=json.load(open(out)) if os.path.exists(out) else []
json.dump(old+res,open(out,'w'))
for r in res: print({k:r.get(k) for k in ('script','env','Iarg','T_0.999','z','rms87','err')})
