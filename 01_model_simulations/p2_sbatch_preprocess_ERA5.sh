#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -p fat_rome
#SBATCH --job-name=p2_sbatch_preproc

# Load modules
module purge
module load 2021

for yr in {2023..2023..1}; do
(
  echo $yr
  conda run -n gtsm3-era5-nrt-slm python p2_preprocess_ERA5.py $yr
) &
done
wait
