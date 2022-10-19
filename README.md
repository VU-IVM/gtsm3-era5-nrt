# gtsm3-era5-nrt

Code for GTSMv3.0 near-real time simulations with ERA5 on Snellius (for Sea Level Monitor)

## About the project
The goal of these scripts are to autmatically run GTSM on Snellius. GTSM is forced with ERA5.

## File explanation

### Bash scripts
The bash scripts are used to define and run the several steps in the workflow. Each step consists of a separate bash script, which calls the relevant python script (or the Delft3D FM singularity container). The workflow.sh contains all steps and does not require specific user input, as the date is automatically retrieved (ensuring all input and output settings are correct).  The following bash scripts called by the workflow.sh script:
workflow --> automate with crontab/job

- ```p1_download_meteo.sh``` takes care of downloading ERA5 and tide data (see p1a download_CDS_ERA5.py and p1b dowload_CDS_tides.py below)
- ```p2_proproc_meteo.sh``` converts the ERA5 data into FM inout format (see p2_preproc_meteo.py below)
- ```p3_prepare_gtsm.sh``` prepares the GTSM run by copying models files from template (see p3_prepare_gtsm_run.py below)
- ```p4_run_gstm.sh``` submits the GTSM run using singularity container
- ```p5_postproc_gtsm.sh``` converts to CDS format and plots the results of the simulations (see p5_postproc_gtsm_run.py below)

### Python scripts
The Python scripts submitted by the bash script are used to download and preprocess the meteorological forcing, to prepare the GTSM simulations, and to postprocess the results. 

- ```p1a download_CDS_ERA5.py``` functionality to download ERA5 using the CDS API (note that a key-file is required)
- ```p1b dowload_CDS_tides.py``` functionality to download GTSM-tides using the CDS API (note that a key-file is required)
- ```p2_preproc_meteo.py``` functionality to convert downloaded data into suitable forcing input for a Delft3D FM model
    - fix lat coordinates and change standard names
    - merge daily files to monthly including spinup
- ```p3_prepare_gtsm_run.py``` functionality to prepare the GTSM model 
    - copy model files and use template to change folder paths 
- ```p5_postproc_gtsm_run.py``` functionality to proprocess GTSM simulation results
    - convert netcdf output into CDS appropiate format 
    - compute residual water levels and annual means
    - plotting results (min, max, mean) on global map
  
To be able to run the Python script you need to install conda and create an virtual environment using the env.yml file (```conda create --name gtsm3-era5-nrt --file env.ym```). This installs the required packages, such as xarray, netCDF4, cartopy, etc. 

## Contact
Sanne Muis - sanne.muis@deltares.nl

More on GTSM: https://publicwiki.deltares.nl/display/GTSM/Global+Tide+and+Surge+Model
More on the Sea Level Monitor: xx
