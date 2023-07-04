#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p thin
#SBATCH -N 1

# load modules
module purge
module load 2021

# loop over months and years
for yr in {1960..1960..1}; do
  for mnth in {1..12..1}; do
  (
    echo $yr $mnth
    conda run -n gtsm3-era5-nrt-slm python p1a_download_ERA5.py $yr $mnth
  ) &
  done
done
wait
