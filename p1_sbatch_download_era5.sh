#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH -N 1 --ntasks-per-node=20
cd /projects/0/GTSM/meteo/ERA5/ecmwf/scripts/
#cp * $TMPDIR
#cd $TMPDIR
module load 2019
module load Miniconda2/4.3.21
source activate ERA5py2env
for i in {1979..1998..1}; do
(
  year="$(printf -- '%04d' "$i" )"
  yearplus="$(printf -- '%04d' "$(($i+1))" )"
  python download_ERA5_CDS.py $year $yearplus
) &
done
wait
for i in {1999..2018..1}; do
(
  year="$(printf -- '%04d' "$i" )"
  yearplus="$(printf -- '%04d' "$(($i+1))" )"
  python download_ERA5_CDS.py $year $yearplus
) &
done
wait
#cp *.nc /projects/0/GTSM/meteo/ERA5/ecmwf/

#!/bin/bash

source=ERA5_hist-cartesius

for yr in {1979..2018..1} #
   do
      name=temp
      rm -f qsub_$name.sh
      echo "#! /bin/sh" >>qsub_$name.sh
      echo "#SBATCH -N 1" >>qsub_$name.sh
      echo "#SBATCH --cpus-per-task=1" >>qsub_$name.sh
      echo "#SBATCH --ntasks=3" >>qsub_$name.sh
      #echo "#SBATCH --ntasks-per-node=1" >>qsub_$name.sh
      echo "#SBATCH --job-name=M_$yr" >>qsub_$name.sh
      echo "#SBATCH -t 35:00" >>qsub_$name.sh
      echo "#SBATCH -o log_$yr.out" >>qsub_$name.sh
      # use python env with netCDF4 and pandas installed
      echo "~/.conda/envs/py3env/bin/python "~/meteo_fields/run_all_spinup.py" $source $yr air_pressure atm msl &" >>qsub_$name.sh
      echo "~/.conda/envs/py3env/bin/python "~/meteo_fields/run_all_spinup.py" $source $yr eastward_wind atm u10 &" >>qsub_$name.sh
      echo "~/.conda/envs/py3env/bin/python "~/meteo_fields/run_all_spinup.py" $source $yr northward_wind atm v10 &" >>qsub_$name.sh
      echo "wait" >>qsub_$name.sh
      chmod +x qsub_$name.sh
      sbatch ./qsub_$name.sh
      rm -f qsub_$name.sh
done



