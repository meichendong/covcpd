#!/bin/sh
script_dir=/pine/scr/m/e/meichen/hovestdta/scripts
cd ${script_dir}
for cc in {1..31}; do
	sed "s|ipath|${cc}|g;"  ${script_dir}/A0009_deng_changepoint_v2.R > sim/a0009/a0009_${cc}.R;
	sbatch -o log/iter/a0009_${cc}.log -n 1 \
		--mem=5000 -t 10:22:00 --job-name="a9_${cc}" \
		/nas/longleaf/apps/r/3.6.0/bin/R --vanilla --no-restore --no-save CMD BATCH \
        	sim/a0009/a0009_${cc}.R sim/a0009/a0009_${cc}.Rout;
done