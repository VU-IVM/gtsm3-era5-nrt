# ---
# Author: N. Aleksandrova
<<<<<<< HEAD
=======
# Contact: natalia.aleksandrova@deltares.nl
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
# Date created: February 2024
# remarks: this script is for collecting all statistical values derived from the GTSM-ERA5-E dataset and assembling it in one netCDF file.

import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import xarray as xr
import numpy as np
import sys
import os
import glob
from datetime import datetime

sys.path.append("..")
from path_dict import path_dict
from pathlib import Path

if __name__ == "__main__":   
  
    # location of EVA and percentiles data
    dir_postproc = path_dict['postproc']
    dir_eva_main = os.path.join(dir_postproc,'EVA-GTSM-ERA5')
    dir_eva = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_1950_2022_1hr_v2') 
    
    #locate .csv files containing extremes
    filenames = 'ds_GTSM-ERA5_1950_2022_stations*eva.csv'
    dir_data = os.path.join(dir_eva,filenames)
    file_list = glob.glob(dir_data)
    file_list.sort()

    #read and merge eva data in one data frame
    ds_gtsm_eva = pd.read_csv(file_list[0])
    for ii in range(1,len(file_list)):
        tmp = pd.read_csv(file_list[ii])
        ds_gtsm_eva = pd.concat([ds_gtsm_eva,tmp])
    del tmp

    #locate and read files that contain coordinates 
    file_list_nc = [str(file_list[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list))]
    ds = xr.open_mfdataset(file_list_nc); ds.close()
    ds = ds.rename({'sea_level_detrended':'wl_quantiles'})
    ds['wl_quantiles'] = ds['wl_quantiles'].astype('float32')

    # add coordinate with return periods
    rps = [1,2,5,10,25,50,75,100]
    ds = ds.assign_coords(return_period=rps)

    # store extreme value analysis MLE model fit parameters 
    ds['eva_param_scale'] = (('stations'), ds_gtsm_eva['scale'].tolist())
    ds['eva_param_shape'] = (('stations'), ds_gtsm_eva['c'].tolist())
    ds['eva_param_loc'] = (('stations'), ds_gtsm_eva['thresh'].tolist())

    # store extreme values
    cols_bf = []; cols_low = []; cols_high = [];
    for rr in rps:
        cols_bf.append(f'{rr}_bf')
        cols_low.append(f'{rr}_lower')
        cols_high.append(f'{rr}_higher')
    ds['wl_extreme_bestfit'] = (('stations','return_period'), ds_gtsm_eva[cols_bf].to_numpy())
    ds['wl_extreme_5perc'] = (('stations','return_period'), ds_gtsm_eva[cols_low].to_numpy())
    ds['wl_extreme_95perc'] = (('stations','return_period'), ds_gtsm_eva[cols_high].to_numpy())

    ds.load()

    ds = ds.rename({'station_x_coordinate':'lon','station_y_coordinate':'lat'})

    # assign coordinate attributes
    ds.station_x_coordinate.attrs = {'units': 'degrees_east', 'standard_name': 'longitude', 'long_name': 'longitude'}
    ds.station_y_coordinate.attrs = {'units': 'degrees_north', 'standard_name': 'latitude', 'long_name': 'latitude'}
    ds.stations.attrs = {'units': '', 'long_name': 'station_id'}
    ds.return_period.attrs = {'units': 'years', 'long_name': 'return_period'}
    ds.wl_quantiles.attrs = {'units': 'm', 'standard_name': 'sea_surface_height_above_mean_sea_level', 'long_name': 'still_water_level_statistics_quantiles', 'description': 'still water level above mean sea level statistical value corresponding to a given quantile based on detrended timeseries over 1950-2022'}
    ds.eva_param_scale.attrs = {'units': '','long_name':'scale_parameter', 'description': 'Scale parameter of the extreme value analysis, POT-GPD method, MLE model'}
    ds.eva_param_shape.attrs = {'units': '','long_name':'shape_parameter', 'description': 'Shape parameter of the extreme value analysis, POT-GPD method, MLE model'}
    ds.eva_param_loc.attrs = {'units': '','long_name':'location_parameter', 'description': 'Location parameter of the extreme value analysis, POT-GPD method, MLE model'}
    ds.wl_extreme_bestfit.attrs = {'units': 'm','standard_name': 'sea_surface_height_above_mean_sea_level','long_name':'extreme_water_level_bestfit', 'description': 'Extreme return value of the still water level based on detrended timeseries over 1950-2022 - best-fit value'}
    ds.wl_extreme_5perc.attrs = {'units': 'm','standard_name': 'sea_surface_height_above_mean_sea_level','long_name':'extreme_water_level_lower_bound', 'description': 'Extreme return value of the still water level based on detrended timeseries over 1950-2022 - lower bound, 5th percentile value'}
    ds.wl_extreme_95perc.attrs = {'units': 'm','standard_name': 'sea_surface_height_above_mean_sea_level', 'long_name':'extreme_water_level_upper_bound', 'description': 'Extreme return value of the still water level based on detrended timeseries over 1950-2022 - higher bound, 95th percentile value'}
    
    # assign global attributes
    ds.attrs={'Conventions':'CF-1.8', 
              'title': 'Statistics of detrended still water level timeseries for GTSM-ERA5-E dataset (1950-2022)', 
              'history': 'This is version 1 of the dataset',
              'institution': 'Deltares', 
              'source': 'GTSMv3 forced with ERA5 climate reanalysis',
              'comment':'',
              'references':'DOI: 10.5281/zenodo.10671284',
              'featureType': 'point', 
              'id': 'GTSM-ERA5-E_water_level_statistics',
              'naming_authority': 'https://deltares.nl/en',
              'summary': 'This dataset has been produced with the Global Tide and Surge Model (GTSM) version 3.0. GTSM was forced with wind speed and pressure fields from ERA5 climate reanalysis covering period 1950-2022. The timeseries at each output point globally were detrended (mean and annual mean removed) and statistics (quantiles and extreme values) were calculated. Extreme values were calculated using POT-GPD method fitted with MLE, using 0.99 percentile as threshold ',
                'date_created': str(datetime.utcnow()) + ' UTC', 
                'date_modified': '', 
                'project': 'Deltares Strategic Research Program', 
                'acknowledgment': 'The development of this dataset was financed with Deltares Strategic Research Program.', 
                'contact': 'natalia.aleksandrova@deltares.nl',
                'license': 'Creative Commons Attribution 4.0 International ', 
                'keywords': 'sea-level rise; climate change; water level; climate; tides; hydrography; global tide and surge model;', 
                'geospatial_lat_min': ds.lat.min().round(3).astype(str).item(), 
                'geospatial_lat_max': ds.lat.max().round(3).astype(str).item(), 
                'geospatial_lon_min': ds.lon.min().round(3).astype(str).item(), 
                'geospatial_lon_max': ds.lon.max().round(3).astype(str).item(), 
                'geospatial_lat_units': 'degrees_north',
                'geospatial_lat_resolution': 'point',
                'geospatial_lon_units': 'degrees_east', 
                'geospatial_lon_resolution': 'point',
                'geospatial_vertical_min': ds['wl_quantiles'].min().round(3).astype(str).item(), 
                'geospatial_vertical_max': ds['wl_extreme_95perc'].max().round(3).astype(str).item(),
                'geospatial_vertical_units': 'm', 
                'geospatial_vertical_positive': 'up',
                'time_coverage_start': '1950-01-01', 
                'time_coverage_end': '2023-12-31'}

    # save dataset
    ofile = Path(os.path.join(dir_eva_main,'ds_GTSM-ERA5-E_1950-2022_stats.nc'))
    ds.to_netcdf(ofile)
