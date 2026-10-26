from bscan import *
if __name__=='__main__':
    jobs=[]
    for key in ('0.045_2','0.07_2'):
        B0=REF[key]['B']
        for I,fs in ((0.0,(0,)),(0.5,(0.004,0.011,0.018)),(1.0,(0.012,0.022,0.032)),(1.5,(0.02,0.032,0.044))):
            for f in fs: jobs.append((key,I,round(B0*(1+f),4)))
    with Pool(4) as p: res=p.map(job,jobs,chunksize=1)
    json.dump(res,open('bscan2.json','w'))
    for r in res: print('%-7s I=%.2f B=%.4f T=%.4f a=%.2f z=%.3f rms=%.2f r85min=%.2f'%(r['key'],r['I'],r['B'],r['T'],r['a']*1e3,r['z'],r['rms87']*1e3,r['r85min']*1e3),flush=True)
