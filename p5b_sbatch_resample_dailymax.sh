#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p fat_rome
#SBATCH --job-name=p5b_sbatch_resample

# load modules
module purge
module load 2021

# settings
scenario=era5
for yr in {1950..1965..1}; do
(
    conda run -n gtsm3-era5-nrt-slm python p5b_resample_TS_DailyMax.py $yr 
) &
done
wait
