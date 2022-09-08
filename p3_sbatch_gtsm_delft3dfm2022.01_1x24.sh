#! /bin/bash
#usage: sbatch ./sbatch_oss68758_1x24.sh
#SBATCH --nodes 1               #-N, -n is total numer of nodes
#SBATCH --ntasks-per-node=24    #you pay for a minimum of 1/4 of the cores on a node
#SBATCH --job-name=gtsm_runname #-J
#SBATCH --time 5-00:00:00       #-t, reduce the expected time if possible to increase your priority
#SBATCH --chdir=./              #chdir set as /path/to/runfolder is useful when calling this script from a different directory
#SBATCH --partition=thin        #type of node
#see sbatch --help for additional parameters
#overview of node types (partition): https://servicedesk.surfsara.nl/wiki/display/WIKI/Snellius+usage+and+accounting
#partition=thin: Single node jobs run on a shared node by default. Add --exclusive if you want to use a node exclusively.
#partition=thin: You will be charged for 0.25 node. You could request 32 CPU cores without increasing the price. Use --ntasks, --cpus-per-task and/or --mincpus.
#partition=thin: Note that by default you will get 1875 MiB of memory per CPU core, unless explicitly overridden by --mem-per-cpu, --mem-per-gpu or --mem.

# stop after an error occured
set -e

#singularity versions from p:\d-hydro\delft3dfm_containers\delft3dfm_2022.01
runscript=/projects/0/einf220/delft3dfm_containers/delft3dfm_2022.01/run_dflowfm_singularity.sh

# load modules
module load 2021
module load intel/2021a

# Partition model, partitioning should be sequential so --ntasks=1, a 2x24 run gives ndomains=48)
$runscript --ntasks 1 --partition:ndomains=24:icgsolver=6 gtsm_model.mdu
# Running model, a 2x24 run gives --nodes=2 (header and srun command), --ntasks=48 (srun command) and --ntasks-per-node=24 (header))
srun --nodes 1 --ntasks 24 $runscript --autostartstop gtsm_model.mdu
