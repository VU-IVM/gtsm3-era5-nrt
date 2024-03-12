#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 128
<<<<<<< HEAD
#SBATCH -p rome
=======
#SBATCH -p fat_rome
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
#SBATCH --job-name=eva_gtsm

# load modules
module purge
module load 2021

step=100
<<<<<<< HEAD
# on large node, 128 cores, it is possible to run 90 batches of 100 stations without running out of memory. When using a different machine, the limits will be different.
# submit analysis of extremes in batches per 100 stations
for st_start in {0..8900..100}; do
(
  st_end=$(( $st_start + $step - 1 ))
  # inputs to the function are start year, end year, station number range and which timeseries to use (1-hour or 10-min, where recommended default is 1-hour due to large memory needed)
  conda run -n gtsm3-era5-nrt-slm python p02_computing_return_periods_GTSM_batch_pyextremes.py 1950 2022 $st_start $st_end "1hr"
=======
# on fat node, 128 cores, it is possible to run 90 batches of 100 stations without running out of memory
for st_start in {0..8900..100}; do
(
  st_end=$(( $st_start + $step - 1 ))
  conda run -n gtsm3-era5-nrt-slm python p02_computing_return_periods_GTSM_batch_pyextremes_serial.py 1979 2018 $st_start $st_end "1hr"
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
) &
done
wait