for k in bscan zs bres mis cor hv; do python3 jobs.py cf_$k.json res_$k.json > log_$k.txt 2>&1; done
