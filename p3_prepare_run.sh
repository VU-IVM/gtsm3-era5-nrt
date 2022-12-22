#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -p thin
#SBATCH --job-name=prepare
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=natalia.aleksandrova@deltares.nl

set -e 

# Load modules
module purge
module load 2021

# Folder paths
base_dir="/gpfs/work1/0/einf3499/"


for yr in {1960..1978..1}; do
(
  echo $yr $base_dir
  conda run -n gtsm-era5-nrt-slm python p3_prepare_run.py $yr $base_dir
) &
done
wait

