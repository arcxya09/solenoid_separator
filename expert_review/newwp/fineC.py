from fine import *
if __name__=='__main__':
    J=[('base',0.03,B,0.3,'src30_0p3.pkl') for B in (2.842,2.856,2.870)]
    J+=[('base',0.03,B,0.5,'src30_0p5.pkl') for B in (2.856,2.870,2.884)]
    J+=[('new',0.06,B,0.75,'src60_0p5.pkl') for B in (1.915,1.925,1.935)]
    J+=[('new',0.06,B,1.0,'src60_1p0.pkl') for B in (1.925,1.932,1.940)]
    with Pool(4) as p:
        for r in p.imap(job,J): print('%-5s I=%.2f B=%.4f T=%.4f a=%.2f z=%.3f rms=%.2f r85min=%.2f %s'%(r['tag'],r['I'],r['B'],r['T'],r['a']*1e3,r['z'],r['rms87']*1e3,r['r85min']*1e3,r['src']),flush=True)
