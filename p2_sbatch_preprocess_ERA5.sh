#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH --tasks-per-node 16

# load modules
module load 2021
conda init bash
conda activate gtsm3-era5-nrt-slm

# settings
pdir="/gpfs/work1/0/einf3499/meteo_ERA5"
spinup=10 # days used for spinup

# yearly runs for 2019, 2020 and 2021
for i in {2019..2021..1}; do
  j=$((i+1))
  (
    tstart=$(date -d $i"-01-01" +'%d-%m-%Y')
    start=$(date -d $i"-01-01" +%Y%m%d)
    sim_start=$(date -d "$start - $spinup days" +%d-%m-%Y)
    sim_end=$(date -d $j"-01-01" +'%d-%m-%Y')
    echo $tstart $sim_start $sim_end
    python p2_preprocess_ERA5_CDS.py $tstart $tend $pdir
  ) &
done
wait