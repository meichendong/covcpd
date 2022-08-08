#!/bin/sh

script_dir=/pine/scr/m/e/meichen/hovestdta/scripts
cd ${script_dir}
for cc in 2; do
	for aa in {1..8}; do					
	sed "s|estorder|${cc}|g; s|ipath|${aa}|g"  ${script_dir}/A0011_cellcycle_cpd.R > sim/a0011/a0011_${aa}_${cc}.R;
			sbatch -o log/iter/a0011_${aa}_${cc}.log -n 1 \
			--mem=8000 -t 08:22:00 --job-name="d_${aa}_${cc}" \
			/nas/longleaf/apps/r/3.6.0/bin/R --vanilla --no-restore --no-save CMD BATCH \
        		sim/a0011/a0011_${aa}_${cc}.R sim/a0011/a0011_${aa}_${cc}.Rout;
	done
done