#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p thin
#SBATCH -N 1

# load modules
module load 2021

<<<<<<< HEAD
pdir="/gpfs/work1/0/einf3499/tides_CDS"
=======
pdir="/gpfs/work1/0/einf3499/CDS_tides"
>>>>>>> 8b09c4ec27384a38290e67c38613e6231798c073
# loop over months and years
for yr in {2018..2025..1}; do
  (
    echo $tstart $pdir
    conda run -n gtsm3-era5-nrt-slm python p1b_download_tides.py $yr $pdir
  ) &
  done
<<<<<<< HEAD
=======
done
>>>>>>> 8b09c4ec27384a38290e67c38613e6231798c073
wait