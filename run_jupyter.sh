#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -p rome

# Make sure the jupyter command is available, either by loading the appropriate modules, sourcing your own virtual environment, etc.
module load 2021

source /home/naleksandro/miniconda3/etc/profile.d/conda.sh
conda activate /home/naleksandro/miniconda3/envs/gtsm3-era5-nrt-slm
 
# Choose random port and print instructions to connect
PORT=`shuf -i 5000-5999 -n 1`
LOGIN_HOST=${SLURM_SUBMIT_HOST}-pub.snellius.surf.nl
BATCH_HOST=$(hostname)
 
echo "To connect to the notebook type the following command from your local terminal:"
echo "ssh -J ${USER}@${LOGIN_HOST} ${USER}@${BATCH_HOST} -L ${PORT}:localhost:${PORT}"
echo
echo "After connection is established in your local browser go to the address:"
echo "http://localhost:${PORT}"

#jupyter notebook --no-browser --port $PORT
jupyter lab --no-browser --port ${PORT}
#jupyter lab --no-browser --port ${port} 