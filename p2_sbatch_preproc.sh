#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -p fat
#SBATCH --job-name=preproc
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=sanne.muis@deltares.nl

set -e 

# Load modules
module purge
module load 2021

# Folder paths
folder_era5="/gpfs/work1/0/einf3499/meteo_ERA5"

# Download ERA5 data and GTSM tides from CDS
conda run -n gtsm3-era5-nrt-slm python p2_preproc_ERA5.py --input_dir $folder_era5 --date_string $1