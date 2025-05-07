#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p rome
#SBATCH --job-name=p5_sbatch_resample

# load modules
module purge
module load 2022

# settings
scenario=era5
for yr in {1950..1950..1}; do
    for mnth in {1..1..1}; do
    (
        conda run -n gtsm3-era5-nrt-slm python p5a_resample_TS.py $yr $mnth
    ) &
    done
done
wait
