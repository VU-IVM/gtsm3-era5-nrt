#!/bin/bash
# Since authentication is required for the SVN repos, this script should be executed and not submitted

# load modules
module purge
module load 2021

#get dir_modelfiles as variable
dir_modelfiles=$(python path_dict.py modelfiles)

mkdir "$dir_modelfiles"
#cd dir_modelfiles

#checkout gtsm3_cmip6 repos folder to dir_modelfiles
svn checkout https://repos.deltares.nl/repos/global_tide_surge_model/trunk/gtsm3_cmip6 "$dir_modelfiles"

#remove non-template ext/mdu/sh files to avoid issues/confusion
rm "$dir_modelfiles"/*.ext
rm "$dir_modelfiles"/*.mdu
rm "$dir_modelfiles"/*.sh
