#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -p thin
#SBATCH --job-name=p4_sbatch_postprocess

# load modules
module purge
module load 2021

# settings
scenario=era5
# yearly runs for 2019, 2020 and 2021
for yr in {1952..1978..1}; do  
(
  conda run -n gtsm3-era5-nrt-slm python 01_postprocess_NL_yearly_stats.py $yr
  conda run -n gtsm3-era5-nrt-slm python 02_postprocess_NL_yearly_TS.py $yr  
) &
done
wait
