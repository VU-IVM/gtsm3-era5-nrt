#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -p rome
#SBATCH -N 1

# load modules
module purge
module load 2021

conda run -n gtsm3-era5-nrt-slm python plot_era5_test.py 

