#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 192
#SBATCH -p genoa
#SBATCH --job-name=eva_gtsm

# load modules
module purge
module load 2022

step=100
# on large node, 128 cores, it is possible to run 90 batches of 100 stations without running out of memory. When using a different machine, the limits will be different.
# submit analysis of extremes in batches per 100 stations
for st_start in {0..1800..100}; do
(
  st_end=$(( $st_start + $step - 1 ))
  # inputs to the function are start year, end year, station number range and which timeseries to use (1-hour or 10-min, where recommended default is 1-hour due to large memory needed)
  conda run -n gtsm3-era5-nrt-slm python p01_computing_return_periods_GTSM_batch_pyextremes.py 1950 2024 $st_start $st_end "1hr"
) &
done
wait