# collector variants x pocket space charge (Rb-87 product current Ip): reflection / capture / landing
# python3 col_sc.py design src B0 out.json
from collector import *
import sys
from multiprocessing import Pool
design,src,B0,out=sys.argv[1],sys.argv[2],float(sys.argv[3]),sys.argv[4]
IPS=[float(x) for x in os.environ.get('IPS','0,1e-4,1e-3').split(',')]
VAR={'pocket300':dict(style='pocket',Vland=300.),
     'pocket1500':dict(style='pocket',Vland=1500.),
     'longdecel1000':dict(style='pocket',Vland=1000.,gap=50e-3),
     'sphere1500':dict(style='sphere',Vland=1500.)}
VARS=os.environ.get('VARS',','.join(VAR)).split(',')
FPK=[float(x) for x in os.environ.get('FPK','0').split(',')]
def job(a):
    vn,Ip,fpk=a; d=DESIGNS[design]; k=1.0 if d['E']<20e3 else 2.5
    gk=dict(VAR[vn]);
    if 'gap' in gk: gk['gap']*=k
    sol=get_sol(design,B0)
    recs,zp,scd=run_beam(design,src,N=int(os.environ.get('NCOL','2000')),B0=B0,I=float(os.environ.get('ITR','0')),iters=4,zplanes=np.arange(d['zc']-0.25,d['zc']+0.01,0.0025))
    col,f7,f5=run_collector(design,recs,zp,sol,Isep=Ip/AB['87'],fpocket=fpk,gk=gk,Efield=[scd['Ef']] if scd.get('Ef') else [],iters=int(os.environ.get('CITERS','5')) if Ip>0 else 1,nmax=int(os.environ.get('NMAX','60000')))
    o=summarize(col,f7,f5,recs['87'].shape[1],recs['85'].shape[1]); o['var']=vn; o['Ip']=Ip; o['fpocket']=fpk; o['P']={kk:(float(v) if not isinstance(v,str) else v) for kk,v in col.P.items()}
    return o
if __name__=='__main__':
    with Pool(int(os.environ.get('NP','4'))) as p: res=p.map(job,[(v,I,f) for v in VARS for I in IPS for f in FPK],chunksize=1)
    json.dump(res,open(out,'w'))
    for o in res:
        a=o['87']; print('fpk=%.2f '%o['fpocket']+'%-14s Ip=%6.1e mA  pocket %.4f refl %.4f cone+dump %.4f  KE %s  ang %s'%(o['var'],o['Ip']*1e3,a.get('pocket',0),a.get('reflected',0),a.get('cone',0)+a.get('dump',0)+a.get('E2',0)+a.get('tube',0)+a.get('shield',0),
              [round(x) for x in a.get('all_KE',[0])],[round(x,1) for x in a.get('all_angle_deg',[0])]))
