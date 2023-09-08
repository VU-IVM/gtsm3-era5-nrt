#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 12
#SBATCH -p rome
#SBATCH --job-name=postprocess_NL

# load modules
module purge
module load 2021

# settings
scenario=era5
for yr in {2022..2022..1}; do  
(
  conda run -n gtsm3-era5-nrt-slm python step1_postprocess_NL_yearly_stats.py $yr
  conda run -n gtsm3-era5-nrt-slm python step2_postprocess_NL_yearly_TS.py $yr  
  conda run -n gtsm3-era5-nrt-slm python step3_postprocess_NL_plots.py $yr  
) &
done
wait
