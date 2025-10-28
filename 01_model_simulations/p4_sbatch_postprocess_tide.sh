#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p rome
#SBATCH --job-name=p4_sbatch_postprocess

# load modules
module purge
module load 2022

# settings - Note: for 2 years of data you need 80 cores on rome partition
scenario=era5
# yearly runs for 2019, 2020 and 2021
for yr in {1950..1950..1}; do
  for mnth in {1..2..1}; do
  (
    conda run -n gtsm3-era5-nrt-slm python p4_postprocess_FM_tide.py $yr $mnth $scenario
  ) &
  done
done
wait
