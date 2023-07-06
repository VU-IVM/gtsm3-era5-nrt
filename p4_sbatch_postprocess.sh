#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p fat
#SBATCH --job-name=p4_sbatch_postprocess

# load modules
module purge
module load 2021

# settings
scenario=era5
# yearly runs for 2019, 2020 and 2021
for yr in {1960..1960..1}; do
  for mnth in {1..12..1}; do
  (
    conda run -n gtsm3-era5-nrt-slm python p4_postprocess_FM.py $yr $mnth $scenario
  ) &
  done
done
wait

# WORKFLOW for postprocessing
#
# 1 remove spinup
# 2 compute surge residual
# 3 change attributes and write monthly files
# 4 ompute and plot min, max, mean for monhtly and annual files 