#!/bin/bash

# Usage:
#   - Modify this script where needed (e.g. number of nodes, number of tasks per node).
#   - Execute this script from the command line of H7 using:
#     sbatch submit_h7.sh
#
# This is an h7 specific script for single or multi-node simulations

#--- Specify Slurm SBATCH directives ------------------------------------------------------------------------
#SBATCH --nodes=1                                 # Number of nodes.
#SBATCH --ntasks-per-node=16                      # The number of tasks to be invoked on each node.
                                                  # For sequential runs, the number of tasks should be '1'.
                                                  # Note: SLURM_NTASKS is equal to "--nodes" multiplied by "--ntasks-per-node".
#SBATCH --job-name=gtsm_era5                      # Specify a name for the job allocation.
#SBATCH --time 0-06:00:00                         # Set a limit on the total run time of the job allocation.
#SBATCH --partition=16vcpu                        # Request a specific partition for the resource allocation.
                                                  # See: https://publicwiki.deltares.nl/display/Deltareken/Compute+nodes.
#SBATCH --mail-type=fail                          # Send an email when the job starts, stops, or fails.
#SBATCH --mail-user=natalia.aleksandrova@deltares.nl      # Specify the email address to which notifications are to be sent.
 
 
# load modules
module load intelmpi/2021.10.0
module load miniforge

# settings
scenario=era5
for yr in {2025..2025..1}; do
    for mnth in {1..12..1}; do
    (
        conda run -n gtsm3-era5 python p5a_resample_TS.py $yr $mnth
    ) &
    done
done
wait
