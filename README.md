# gtsm3-era5-nrt

Code for GTSMv3.0 near-real time simulations with ERA5 on Snellius (for Sea Level Monitor)

## About the project
The goal of these scripts are to autmatically run GTSM on Snellius. GTSM is forced with ERA5.

## File explanation

### Bash scripts
The bash scripts are used to define and run the several steps in the workflow. Each step consists of a separate bash script, which calls the relevant python script (or the Delft3D FM singularity container). The workflow.sh contains all steps and does not require specific user input, as the date is automatically retrieved (ensuring all input and output settings are correct).  The following bash scripts called by the workflow.sh script:
workflow --> automate with crontab/job

- '''p1_download_meteo.sh takes care of downloading ERA5 and tide data
- '''p2_proproc_meteo.sh''' converts the ERA5 data into FM inout format (see download_data.py and convert_data.py below)
- '''p3_prepare_gtsm.sh''' prepares the GTSM run by copying models files from template
- '''p4_run_gstm.sh''' submits the GTSM run using singularity container
- ''''p5_postproc_gtsm.sh''' converts to CDS format and plots the results of the simulations    

### Python scripts

    download_data.py functionality to download ERA5 and SEAS5 data using the CDS API (note that a key-file is required)
    convert_data.py functionality to convert downloaded data into suitable forcing input for a wflow_sbm model
    plot_wflow_results.py functionality to convert simulation results into a figure showing the last month of discharge, together with the forecasted discharge.

## Usage

## Contact

Sanne Muis - sanne.muis@deltares.nl

More on GTSM: https://publicwiki.deltares.nl/display/GTSM/Global+Tide+and+Surge+Model
More on the Sea Level Monitor: xx




# Step 1: Download data
- 1a download daily ERA5 files from CDS
- 1b dowload yearly tide files from CDS

# Step 2: Preprocess ERA5 to FM format
- fix lat coordinates
- change standard names
- merge daily files to yearly/monthly

# Step 3: Start GTSM model simulation
- use singularity containers
- template files

# Step 4: Postprocessing of raw data
- compute residual water levels
- fix attributes
- compute annual means
- plot data
