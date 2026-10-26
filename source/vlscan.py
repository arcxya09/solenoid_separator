from scrun import *
from multiprocessing import Pool
def j(v): 
    r=run(0,v,N=800,seed=3); return r
with Pool(4) as p: res=p.map(j,[-12000,-13000,-14000,-15000,-16000,-17000,-18000,-19000])
json.dump(res,open('vlscan.json','w'))
for r in res: print(r['Vl'],round(r['rp_rms']*1e3,2),round(r['rp_mean']*1e3,2),round(r['rwid_rms']*1e3,2))
