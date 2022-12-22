#!/bin/bash
#SBATCH -t 002:00:00
#SBATCH -N 1
#SBATCH -p thin
#SBATCH --job-name=prepare
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=sanne.muis@deltares.nl

set -e 

# Load modules
module purge
module load 2021

# Folder paths
base_path="/gpfs/work1/0/einf3499/"

# Download ERA5 data and GTSM tides from CDS
conda run -n gtsm3-era5-nrt-slm python p3_prepare_run.py --base_dir $base_path --date_string $1