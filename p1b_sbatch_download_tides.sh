#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p thin
#SBATCH -N 1

# load modules
module load 2021

pdir="/gpfs/work1/0/einf3499/tides_CDS"
# loop over months and years
for yr in {2018..2025..1}; do
  (
    echo $tstart $pdir
    conda run -n gtsm3-era5-nrt-slm python p1b_download_tides.py $yr $pdir
  ) &
  done
wait