from hc30 import *
if __name__=='__main__':
    J=[(B,I,'src104_30k_5mA.pkl',1500) for I,Bs in ((6,(1.928,1.936)),(7,(1.932,1.940)),(8,(1.935,1.945))) for B in Bs]
    with Pool(4) as p:
        for r in p.imap(job,J): print('I=%g B=%.4f T=%.4f a=%.2f z=%.3f rms=%.2f r85min=%.2f'%(r['I'],r['B'],r['T'],r['a']*1e3,r['z'],r['rms87']*1e3,r['r85min']*1e3),flush=True)
