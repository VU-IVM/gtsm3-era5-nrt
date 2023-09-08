#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p thin
#SBATCH --job-name=eva_gtsm

# load modules
module purge
module load 2021

conda run -n gtsm3-era5-nrt-slm python p02_computing_return_periods_GTSM.py 1985 2014

