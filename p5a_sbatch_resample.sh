#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p fat
#SBATCH --job-name=p5_sbatch_resample

# load modules
module purge
module load 2021

# settings
scenario=era5
for yr in {1951..1960..1}; do
    for mnth in {1..12..1}; do
    (
        conda run -n gtsm3-era5-nrt-slm python p5a_resample_TS.py $yr $mnth
    ) &
    done
done
wait
