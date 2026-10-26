import tol
from tol import *
tol.SIG.update(src_tilt=1.0e-3,b2=2.5e-4)
if __name__=='__main__':
    nom=json.load(open('tol1.json'))['nom']; A={c:v['a'] for c,v in nom.items()}
    J=[(c,i,1.0,A[c]) for c in ('new0','new05') for i in range(20)]
    with Pool(4) as p: mc=p.map(mcjob,J,chunksize=1)
    json.dump(dict(SIG=tol.SIG,mc=mc),open('tol2.json','w'))
    for c in ('new0','new05'):
        for lev in ('L0','L1'):
            T=np.array([r[lev]['T'] for r in mc if r['case']==c]); P=np.array([r[lev].get('pur',1.0) for r in mc if r['case']==c])
            print(c,lev,'mean T %.4f P10 %.4f min %.4f minpur %.5f yield %.2f'%(T.mean(),np.percentile(T,10),T.min(),P.min(),((T>=0.99)&(P>=0.999)).mean()),flush=True)
