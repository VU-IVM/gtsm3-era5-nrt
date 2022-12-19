#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p thin
#SBATCH -N 1

# load modules
module purge
module load 2021

pdir="/gpfs/work1/0/einf3499/tides_CDS_extended"

# loop over months and years
for i in {1952..1978..1}; do
  for j in {1..12..1}; do
  (
    tstart=$(date -d ""$(printf -- '%04d' "$i" )"-"$(printf -- '%02d' "$j" )"-01" +'%d-%m-%Y')
    echo $tstart $pdir
    conda run -n gtsm-era5-nrt-slm python p1b_download_tides.py $tstart $pdir
  ) &
  done
done
wait

