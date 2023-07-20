#! /bin/bash
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --job-name=%JOBNAME%
#SBATCH --time 5-00:00:00

# purpose of run
export purpose="GTSMv3.0 - ERA5 run near-realtime for Sea Level Monitor"
echo "========================================================================="
echo "Submitting Dflow-FM run $name in  $PWD"
echo "Purpose: $purpose"
echo "Starting on $SLURM_NTASKS domains, SLURM_NNODES nodes"
echo "Wall-clock-limit set to $maxwalltime"
echo "========================================================================="

# stop after an error occured
set -e

# load modules
module purge
module load 2021
module load intel/2021a

#singularity versions
modelFolder=${PWD}/../../..
singularityFolder=/projects/0/einf3499/99_delft3dfm_containers/delft3dfm_2022.04/

# DIMR input-file; must already exist!
dimrFile=dimr_config.xml
mduFile=gtsm_fine.mdu

# Replace number of processes in DIMR file
PROCESSSTR="$(seq -s " " 0 $((SLURM_NTASKS-1)))"
sed -i "s/\(<process.*>\)[^<>]*\(<\/process.*\)/\1$PROCESSSTR\2/" $dimrFile

# Replace MDU file in DIMR-file
sed -i "s/\(<inputFile.*>\)[^<>]*\(<\/inputFile.*\)/\1$mduFile\2/" $dimrFile

# Partition model, partitioning should be sequential so --ntasks=1, a 2x24 run gives ndomains=48
$singularityFolder/execute_singularity_snellius.sh $modelFolder run_dflowfm.sh --partition:ndomains=$SLURM_NTASKS:icgsolver=6 $mduFile
# Running model, a 2x24 run gives --nodes=2, --ntasks=48 and --ntasks-per-node=24
$singularityFolder/execute_singularity_snellius.sh $modelFolder run_dimr.sh -m $dimrFile