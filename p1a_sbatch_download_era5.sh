#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p thin
#SBATCH -N 1

# load modules
module load 2021

pdir="/gpfs/work1/0/einf3499/meteo_ERA5_extended"
# loop over months and years
for yr in {1979..1979..1}; do
  for mnth in {1..12..1}; do
  (
    echo $yr $mnth $pdir
    conda run -n gtsm3-era5-nrt-slm python p1a_download_ERA5.py $yr $mnth $pdir
  ) &
  done
done
wait