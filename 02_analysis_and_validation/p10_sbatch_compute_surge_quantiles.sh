#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p rome
#SBATCH --job-name=surge_q

# load modules
module purge
module load 2022

conda run -n gtsm3-era5-nrt-slm python p10_compute_surge_quantiles.py 
