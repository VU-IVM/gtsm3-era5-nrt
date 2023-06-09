#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -p thin
#SBATCH --job-name=p3_prepare_run

# Load modules
module purge
module load 2021

for yr in {1960..1978..1}; do
(
  echo $yr
  conda run -n gtsm3-era5-nrt-slm python p3_prepare_run.py $yr
) &
done
wait
