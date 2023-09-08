#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -p fat
#SBATCH --job-name=eva_gtsm

# load modules
module purge
module load 2021

conda run -n gtsm3-era5-nrt-slm python p02_computing_return_periods_GTSM.py 1951 2022

conda run -n gtsm3-era5-nrt-slm python p02_computing_return_periods_GTSM.py 1979 2018