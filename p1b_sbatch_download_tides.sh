#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p thin
#SBATCH -N 1

# load modules
module purge
module load 2021

# loop over months and years
for yy in {1950..1950..1}; do
  for mm in {2..5..1}; do
  (
    yr=$(printf '%i' "$yy")
    mnth=$(printf '%i' "$mm")
    echo $yr $mnth
    conda run -n gtsm3-era5-nrt-slm python p1b_download_tides.py $yr $mnth
  ) &
  done
done
wait
