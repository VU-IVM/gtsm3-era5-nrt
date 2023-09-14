#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p fat
#SBATCH --job-name=eva_gtsm

# load modules
module purge
module load 2021

step=1000

for st_start in {0..11000..1000}; do
(
  st_end=$(( $st_start + $step - 1 ))
  conda run -n gtsm3-era5-nrt-slm python p02_computing_return_periods_GTSM_batch.py 1951 2022 $st_start $st_end "1hr"
) &
done
wait