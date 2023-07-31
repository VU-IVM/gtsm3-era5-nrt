#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 12
#SBATCH -p thin
#SBATCH --job-name=postprocess_NL

# load modules
module purge
module load 2021

# settings
scenario=era5
# yearly runs for 2019, 2020 and 2021
for yr in {1952..1978..1}; do  
(
  conda run -n gtsm3-era5-nrt-slm python step1_postprocess_NL_yearly_stats.py $yr
  conda run -n gtsm3-era5-nrt-slm python step2_postprocess_NL_yearly_TS.py $yr  
) &
done
wait
