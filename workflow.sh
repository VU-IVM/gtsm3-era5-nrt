# Runs the gtsm workflow as one script, based on last month

scriptsdir="/gpfs/work1/0/einf3499/scripts"

d=$(date)
thisym=$(date "+%Y_%m")
lastym=$(date --date "$d -1 month" "+%Y_%m") 

## use set -e in bash script to catch error
## use crontab / scrontab

# 1 Download daily ERA5 data and tides 
downloadjobid=$(sbatch $scriptsdir/p1a_sbatch_download_era5.sh $lastym | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started download job with id: $downloadjobid"
# 2 Prepare data and convert to FM format
preprocjobid=$(sbatch --dependency=afterok:$downloadjobid $scriptsdir/p2_sbatch_preprocess_ERA5.sh $lastym | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started preproc job with id: $preprocjobid"
# 3 Prepare run 
preparejobid=$(sbatch --dependency=afterok:$preprocjobid $scriptsdir/p3_prepare_and_submit_runs_ERA5_yearly.py$lastym | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started prepare job with id: $preparejobid"
# 4 Run GTSM --> check folder
runjobid=$(sbatch --dependency=afterok:$preparejobid $scriptsdir/p3_sbatch_gtsm_delft3dfm2022.01_xnodes.sh $lastym | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started run job with id: $runjobid"
# 5 Postprocess and plot
postprocjobid=$(sbatch --dependency=afterok:$catchupjobid $scriptsdir/p4_sbatch_postprocess.sh $lastym | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started postproc job with id: $postprocjobid"
