from fine import *
if __name__=='__main__':
    J=[('newV3',0.06,round(1.894*(1+f),4),I,'src_v3.pkl') for I in (1.0,1.5) for f in (0.01,0.02,0.03)]
    with Pool(4) as p:
        for r in p.imap(job,J): print('%-5s I=%.2f B=%.4f T=%.4f a=%.2f z=%.3f rms=%.2f r85min=%.2f'%(r['tag'],r['I'],r['B'],r['T'],r['a']*1e3,r['z'],r['rms87']*1e3,r['r85min']*1e3),flush=True)
