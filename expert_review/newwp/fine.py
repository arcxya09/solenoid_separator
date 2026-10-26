import wp
from wp import *
from bscan import stage_any
from multiprocessing import Pool
def job(a):
    tag,R0,B,I,src=a; wp.SRC=src; nm='f%s_%g_%g'%(tag,B,I); mk(nm,R0,B,'C3')
    recs,zp,info=run(nm,N=1500,I=I*1e-3,zplanes=np.arange(1.25,1.5,0.005),seed=1,iters=4)
    s=stage_any(recs,zp); s.update(tag=tag,I=I,B=B,src=src); return s
if __name__=='__main__':
    J=[]
    S={0:'src60_0p2.pkl',0.5:'src60_0p5.pkl',1.0:'src60_1p0.pkl',1.5:'src60_1p5.pkl'}
    for I,fs in ((0,(-0.01,0,0.01)),(0.5,(0,0.01,0.02,0.03)),(1.0,(0,0.01,0.02,0.03)),(1.5,(0.01,0.02,0.03,0.04))):
        for f in fs: J.append(('new',0.06,round(1.894*(1+f),4),I,S[I]))
    for I,fs in ((0,(-0.01,0,0.01)),(0.2,(0,0.01,0.02,0.03)),(0.3,(0.01,0.02,0.03,0.04)),(0.5,(0.02,0.03,0.04,0.05))):
        for f in fs: J.append(('base',0.03,round(2.80*(1+f),4),I,'src_v3.pkl'))
    with Pool(4) as p:
        res=[]
        for r in p.imap(job,J):
            res.append(r); print('%-5s I=%.2f B=%.4f T=%.4f a=%.2f z=%.3f rms=%.2f r85min=%.2f'%(r['tag'],r['I'],r['B'],r['T'],r['a']*1e3,r['z'],r['rms87']*1e3,r['r85min']*1e3),flush=True)
    json.dump(res,open('fine.json','w'))
