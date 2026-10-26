python3 -c "
from bl import *
s=get_sol(2.785); out=run(s,N=10000); np.savez('base2785.npz',o85=out['85'],o87=out['87']); print(best(out))" > log_base2785.txt 2>&1
for k in mis cor hv bres; do python3 jobs.py cf2_$k.json res2_$k.json > log2_$k.txt 2>&1; done
python3 sc3.py cf_sc.json res_sc.json env > log_sc.txt 2>&1
python3 scat.py > log_scat.txt 2>&1
