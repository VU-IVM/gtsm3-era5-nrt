#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -p fat
#SBATCH --job-name=p5_sbatch_resample

# load modules
module purge
module load 2021

# settings
scenario=era5
for yr in {2019..2022..1}; do
(
    conda run -n gtsm3-era5-nrt-slm python p5b_resample_TS_yearly.py $yr
) &
done
wait
