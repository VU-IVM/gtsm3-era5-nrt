#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=sanne.muis@deltares.nl

set -e 

# Load modules
module purge
module load 2021

# Folder paths
folder_era5="/gpfs/work1/0/einf3499/meteo_ERA5"
folder_tides="/gpfs/work1/0/einf3499/tides_CDS"

# Download ERA5 data and GTSM tides from CDS
conda run -n gtsm3-era5-nrt-slm python p1a_download_ERA5.py --output_dir $folder_era5 --date_string $1
conda run -n gtsm3-era5-nrt-slm python p1b_download_tides.py --output_dir $folder_tides --date_string $1