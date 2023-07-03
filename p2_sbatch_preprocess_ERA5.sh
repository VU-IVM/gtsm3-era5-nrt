#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH -p fat
#SBATCH --job-name=p2_sbatch_preproc

# Load modules
module purge
module load 2021

for yr in {1950..1951..1}; do
(
  echo $yr
  conda run -n gtsm-era5-nrt-slm python p2_preprocess_ERA5.py $yr
) &
done
wait

# WORKFLOW for preparation of forcing files
#
# yearly runs (or shorter) are mandatory due to the current implementation of sea-level rise, where the model is initialised with the given SLR and not updated (FINITE VOLUME) over time.
#
# 1. METEO-FORCING (yearly):
#	-Modify for correct longitude range [-180 to 180] and overlap
#	-right var names and attributes
#	-include spinup in yearly files
#	-Set initial timesteps to zero to allow for SLR correction
# 2. SLR --> files are provided in the folder GTSMip\forcing_files\slr
#	-Modify for correct longitude range [-180 to 180] and overlap
#	-translate to pressure
#	-right var names and attributes
#	-include spinup in yearly files - this we didn't do like this before, but constant yearly files would be more robust, now the pressure field will vary in time.
#	- should be included in .ext file
# 3.Mean sea level --> files are provided in the folder GTSMip\forcing_files\meteo
#	- this files provides a correction of the mean sea level based on the mean sea level pressure over 1986-2005 in ERA-Interim, this to reference the output to MSL
#	- should be included in .ext file