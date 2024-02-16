#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -p fat_rome
#SBATCH --job-name=eva_gtsm

# load modules
module purge
module load 2021

conda run -n gtsm3-era5-nrt-slm python p03_plot_eva_GTSM_comparison.py 
