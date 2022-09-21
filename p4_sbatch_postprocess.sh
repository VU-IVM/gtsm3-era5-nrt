#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -p fat
#SBATCH --job-name=preprocERA5

# load modules
module load 2021

# settings
pdir=/gpfs/work1/0/einf3499/slr_tide_surge_runs/
scenario=era5
# yearly runs for 2019, 2020 and 2021
for yr in {2019..2021..1}; do
  ( sim_start=$(date -d $i"-01-01" +%Y%m%d)
    echo $sim_start 
    conda run -n gtsm3-era5-nrt-slm python p4_postprocess_FM.py $yr $senario $pdir
    ) &
done
wait

# WORKFLOW for postprocessing
#
# 1 remove spinup
# 2 compute surge residual
# 3 change attributes and write monthly files
# 4 ompute and plot min, max, mean for monhtly and annual files 