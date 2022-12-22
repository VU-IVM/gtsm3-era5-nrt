#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -p fat
#SBATCH --job-name=preproc
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=natalia.aleksandrova@deltares.nl

set -e 

# Load modules
module purge
module load 2021

# Folder paths
folder_era5="/gpfs/work1/0/einf3499/meteo_ERA5_extended"

for yr in {1952..1969..1}; do
(
  echo $yr $folder_era5
  conda run -n gtsm-era5-nrt-slm python p2_preproc_ERA5.py $yr $folder_era5
) &
done
wait

