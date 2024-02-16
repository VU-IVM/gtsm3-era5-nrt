#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p fat_rome
#SBATCH --job-name=eva_gtsm

# load modules
module purge
module load 2021

step=100
# on fat node, 128 cores, it is possible to run 90 batches of 100 stations without running out of memory
for st_start in {0..8900..100}; do
(
  st_end=$(( $st_start + $step - 1 ))
  conda run -n gtsm3-era5-nrt-slm python p02_computing_return_periods_GTSM_batch_pyextremes_serial.py 1979 2018 $st_start $st_end "1hr"
) &
done
wait