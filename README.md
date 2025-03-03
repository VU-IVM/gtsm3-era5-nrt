# GTSM-ERA5-E: an extended reanalysis dataset of global water levels

## About the project
Extreme sea levels, generated by storm surges and high tides and increasing with sea level rise, have the potential to cause coastal flooding and erosion. Global datasets are instrumental for mapping extreme sea levels and associated societal risks. 

Harnessing the backward extension of the [ERA5 reanalysis](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview), we produced a dataset containing timeseries of tides and storm surges based on a global hydrodynamic model that covers the period **1950-2024**. The dataset is an extension of a previously published dataset that covered a shorter period (1979-2018). Using the extended dataset, we estimated extreme sea levels for various return periods globally. Validation showed that there is a good agreement between observed and modelled extreme sea levels. The extended 75-year dataset allows for a more robust estimation of return periods, often resulting in smaller uncertainties than its 40-year precursor. This underscores the necessity for long timeseries and highlights the strength of long-term modelling enabled by the ERA5 reanalysis database extension. The present dataset can be used for assessing flood risk, climate variability and climate changes. 

The background of the dataset will be described in more detail in the paper by Aleksandrova et al. (2024) [^1] (currently under review).

## About this repository
The goal of the project was to develop a reproducible workflow used to generate datasets of water level and surge timeseries, and to compute extreme water levels globally. The scripts in this repository were used to produce the GTSM-ERA5-E dataset and include the data processing necessary for the hydrodynamic model simulations, subsequent statistical data analysis and validation against observations. The different steps of the overall workflow are shown below.

<img src="images/GTSM-ERA5-workflow.png" width=60%>

### Preparation
To be able to run the Python scripts you need to install conda and create an virtual environment using the environment.yml file (`conda env create --file environment.yml`). This installs the required packages, such as xarray, netCDF4, cartopy, etc. You also need to specify the relevant folders in the path_dict.py file. 

## Model simulations
The main part of the workflow are the model simulations, which are used to produce hourly timeseries of total water levels and surge levels globally. We make use of the same modelling approach as presented in Muis et al. (2020,2024)[^2][^3]. We further develop a semi-automated and portable workflow that can easily be used to deploy Global Tide and Surge Model (GTSM) on a high-performance computing cluster. In this case we made use of the [Dutch National Supercomputer Snellius](https://www.surf.nl/en/dutch-national-supercomputer-snellius) by SURF. Snellius makes use of SLURM for managing and scheduling jobs on Linux clusters. 

### Hydrodynamic modelling software
The modelling of global tides and surges is done by using the GTSMv3.0 model. The [Global Tide and Surge Model (GTSM)](https://publicwiki.deltares.nl/display/GTSM/Global+Tide+and+Surge+Model) is a depth-averaged hydrodynamic model, developed by [Deltares](https://www.deltares.nl/), with global coverage. GTSM is based on [Delft3D Flexible Mesh software](https://www.deltares.nl/en/software-and-data/products/delft3d-fm-suite/modules/d-flow-flexible-mesh) and has a spatially varying resolution which increases towards the coast. GTSM can be used to simulate water levels and currents, that arise from tides and storm surges. The model has showed to be able to simulate tides and storm surge with enough accuracy when forced with wind and pressure fields from the ERA5 climate reanalysis[^4][^5].  

To run the model, we make use of a Delft3D FM Singularity container, which allows us to run Delft3D FM on any HPC clusters in a simple, portable, and reproducible way. More instructructions on how to obtain a delft3dfmcontainer can be found [here](https://oss.deltares.nl/web/delft3dfm/get-started). Deltares collegues can go to `p:\d-hydro\delft3dfm_containers` to obtain a container. The `fm_container` folder can be adjusted in `path_dict.py`.

### Bash scripts
The model simulation consists of a sequence of bash scripts with SLURM commands that calls the relevant python scripts (or the Delft3D FM singularity container in case of model simulations).

- `p1a_sbatch_download_era5.sh` - for downloading ERA5 data (see p1a download_CDS_ERA5.py below)
- `p1b_sbatch_download_tides.sh` - for downloading tidal data (see p1b dowload_CDS_tides.py below)
- `p1c_checkout_gtsm3_cmip6.sh` - checkout gtsm3_cmip6 repos to the `modelfiles` folder, this contains the files needed to run the GTSM model and template files for model settings 
- `p2_sbatch_preprocess_ERA5.sh` - converts the ERA5 data into FM input format (see p2_sbatch_preprocess_ERA5.py below)
- `p3_prepare_run.sh` - prepares the GTSM run by copying model input files and adjusting the template (see p3_prepare_run.py below)
- `p4_sbatch_postprocess.sh` - converts model output to timeseries files containing 10-min total water levels and surge levels (see p4_postprocess_FM.py below)
- `p5a_sbatch_resample.sh` - resamples 10-min timeseries to hourly timeseries (see p5a_resample_TS.py below)
- `p5b_sbatch_resample_dailymax.sh` - resamples 10-min timeseries of total water levels to daily maxima timeseries (see p5b_resample_TS_DailyMax.py below)

### Python scripts
The Python scripts submitted by the bash script are used to download and preprocess the meteorological forcing, to prepare the GTSM simulations, and to postprocess the results. 

- `p1a download_ERA5.py` - functionality to download ERA5 data (atm. pressure and wind fields) using the CDS API (note that a key-file is required)
- `p1b dowload_tides.py` - functionality to download GTSM-derived tidal level timeseries using the CDS API (note that a key-file is required)
- `p2_preprocess_ERA5.py` - functionality to convert downloaded data into correct format of forcing input for the GTSM (Delft3D-FM) model. This includes prescribing correct longitude ranges, merging yearly forcing data in one file with spinup period included, ensuring correct variable names and attributes.
- `p3_prepare_run.py` - functionality to prepare the GTSM model runs by copying model input files and adjusting template files with model settings in the directory with model runs
- `p4_postprocess_FM.py` - functionality to postprocess the data from the GTSM model results into 10-min timeseries data files for total water levels and surge levels with appropriate attributes
- `p5a_resample_TS.py` - functionality to resample 10-min timeseries data to hourly timeseries (total water level and surge).
- `p5b_resample_TS_DailyMax.py` - functionality to resample 10-min timeseries data to daily maxima timeseries (total water levels only).

## Statistical analysis and validation
Based on model output results, statistical analysis is performed on the timeseries from 1950 to 2022. This includes computing the various percentiles and performing extreme value analysis for all output points of the GTSM model (43000+ stations). Validation of the timeseries is performed by comparing it to the Global Extreme Sea Level Analysis (GESLA) dataset[^6]. A combination of a bash and Python script is used to download the GTSM-ERA5 water level global reanalysis timeseries for years 1979-2018, and to perform extreme value analysis on the full timeseries. The GTSM-ERA5 which was already available on [CDS](https://cds.climate.copernicus.eu/cdsapp#!/dataset/sis-water-level-change-timeseries-cmip6?tab=overview) - this dataset was generated using the same methodology and is extended in the current project. The scripts used for plotting and validation were ran without bash scripts, but using an interactive computational node with Jupyter Lab. 

### Scripts
- `p00_sbatch_download_GTSM.sh` and `p00_download_GTSM_waterlevels_from_CDS.py`- scripts to download the GTSM-ERA5 water level timeseries dataset for 1979-2018 from GTSM
- `p01_sbatch_compute_EVA_GTSM.sh` and `p01_computing_return_periods_GTSM_batch_pyextremes.py` - scripts to compute statistics (percentiles and extreme values) based on water level timeseries for a chosen period (we did this analysis on 1950-2022 and on 1979-2018 timeseries)
- `p02_plot_eva_GTSM_comparison.py` - script used for plotting percentiles and extreme values of total water levels on a global map, comparison between the GTSM-ERA5-E (1950-2022) and GTSM-ERA5 (1979-2018) datasets, and plots of timeseries and statistics for a few chosen locations around the globe.
- `p03_plot_eva_GTSM_GESLA.py` - script used for comparing model data to the observational data from GESLA dataset, and making respective plots. This script uses a subset of the GESLA dataset as input.
- `p04_make_dataset_GTSM-ERA5-E_stats.py` - script used to generate a dataset file containing the statistics (percentiles and extreme values) of total water levels based on 1950-2022 timeseries for all GTSM output locations.

## Data sources
We make use of various data sources. Whereas most of them are open and retrieval is automated and part of the scripts, there are some data sources for which retrieval is not automated yet. This because there are not available in open repositories or have been produced specifically for this project. These sources are described below. 

### ERA5
We use the ERA5 climate renalysis  as meteorological forcing (`10m_u_component_of_wind`,`10m_v_component_of_wind` and `mean_sea_level_pressure`). The download of this data using [Copernicus Climate Data Store API](https://doi.org/10.24381/cds.adbb2d47) is part of the worflow. 

### Tides
Tides (`tidal_elevation`) are downloaded from the [Copernicus Climate Data Store API](https://doi.org/10.24381/cds.a6d42d60). The tides are computed using the same modelling approach based on GTSMv3.0 as the total water level simulations and are used to compute the storm surge levels. 

### Mean sea level and vertical reference 
The mean sea level is prescribed in yearly increments to include sea level rise. The mean sea level files have been produced by Dewi Le Bars from KNMI. The period 1986-2005 is used as reference period. Sea level fields are computed from the sum of different contributors, including dynamic changes, thermal expansion, changes in gravitational fields, and contribution from glaciers and ice sheets. The different contributions are computed and combined using the probabilistic model described in Le Bars (2018). For the period 1950-2016, we use products based on observations for the Antarctic and Greenland ice sheets (Mouginot et al., 2019; Rignot et al., 2019), the glaciers (Marzeion et al., 2015), thermal expansion between 0 and 2000 m depth (Levitus et al., 2012), and climate-driven water storage (Humphrey & Gudmundsson, 2019). The ice sheets are assumed to be in equilibrium before 1979 for Antarctica and 1972 for Greenland because no data are available before these dates. For the period 2016-2050 we use sea-level rise projections based on the Fifth Assessment Report (AR5) of the Intergovernmental Panel on Climate Change (IPCC) for the RCP8.5 scenario (Church et al., 2013), very similar to the SSP5-8.5 scenario used by the models as above. The redistribution of water in the ocean due to wind changes and local steric effects is taken from the CMIP5 models (i.e. ‘zos’ field for the entire period). The fingerprints for the ice sheets, glaciers and land water storage are from the AR5 assessment, and include the gravitational, rotational and Earth elastic response. For the dynamics of the Antarctic contribution we use the re-evaluation presented in the IPCC’s Special Report on the Ocean and Cryosphere in a Changing Climate (SROCC) (Oppenheimer et al., 2019). Additionally, we add the glacial isostatic adjustment from the ICE-6G model (Peltier et al., 2015) but do not consider other processes of vertical land motion, such as subsidence or tectonics. The uncertainty in mean sea level is removed by selecting the median of the sea level observations and projections distributions. 

The mean sea level rise (relative to 1986-2005) fields were converted to a pressure field that is used as an input to the GTSM model runs. For simulations with atmospheric forcing (ERA5), a correction on pressure fields is introduced by subtracting the mean sea-level pressure field (MSLP) over 1986–2005 (based on ERA-Interim because ERA5 was not available at the time when these input files were created). This is done to maintain the MSL in the model as relative to 1986-2005 mean. Essentially, the pressure field in the GTSM-ERA5 simulations consists of three components that are summed up: ERA5 atmospheric pressure, MSL correction (1986-2005) and pressure correction corresponding to the sea level rise. At model initialization, the first two components are set to 0, which results in the correct initial mean sea level in the model simulation for the given year (correct total volume of water). This approach is designed to mimic sea level rise within the GTSM model.

The bathymetry dataset in the GTSMv3.0 model is based on GEBCO_2014 (Muis et al., 2020)[^2], which is referenced to MSL. While there are newer (and possibly more accurate) global bathymetry datasets from GEBCO, we choose to use this dataset to stay consistent with the previous GTSM-ERA5 runs and to keep the same bathymetry source that the GTSMv3.0 model was calibrated with.

Both mean sea level and vertical reference files can be downloaded [here](https://doi.org/10.5281/zenodo.3948088).

### Validation data
Validation is performed using the observational data from the [GESLA (Global Extreme Sea Level Analysis) project](https://www.gesla.org), version 3 of the dataset. A subset of this dataset was made to include only stations where the available observations span a period of at least 50 years overlapping with the period covered by GTSM-ERA5-E dataset (1950-2022), and where no more than 25% of data is missing from that period. Within the script used for validation an additional filter is applied to only use stations located no further than 10 km from the nearest model output locations.

### Previous GTSM-ERA5 dataset
The previous GTSM-ERA5 reanalysis dataset, covering 1979-2018, is available at the [C3S Climate Data Store](https://doi.org/10.24381/cds.a6d42d60).  

## Contact
Sanne Muis - sanne.muis@deltares.nl  
Natalia Aleksandrova - natalia.aleksandrova@deltares.nl

## References
[^1]: Aleksandrova, N., Veenstra, J., Gwee, R. & Muis, S. (2024). Global dataset of storm surges and extreme sea levels for 1950-2022 based on the ERA5 climate reanalysis. In review.
[^2]: Muis, S., Apecechea, M. I., Dullaart, J., ... & Verlaan, M. (2020). A High-resolution global dataset of extreme sea levels, tides, and storm surges, including future projections. Frontiers in Marine Science, doi:10.3389/fmars.2020.00263
[^3]: Muis, S. et al (2023). Global projections of storm surges using high-resolution CMIP6 climate models. In review.
[^4]: Muis, S., Verlaan, M., Winsemius, H. C., Aerts, J. C. J. H., & Ward, P. J. (2016). A global reanalysis of storm surge and extreme sea levels. Nature Communications, doi:10.1038/ncomms11969.
[^5]: Dullaart, J. C., Muis, S., Bloemendaal, N., & Aerts, J. C. (2020). Advancing global storm surge modelling using the new ERA5 climate reanalysis. Climate Dynamics, doi:10.1007/s00382-019-05044-0
[^6]: Haigh, I.D. et al. GESLA Version 3: A major update to the global higher-frequency sea-level dataset. Geoscience Data Journal (2022) doi:10.1002/gdj3.174.


