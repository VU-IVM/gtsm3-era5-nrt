# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 2023
@author: n-aleksandrova,veenstra

This script is for resampling the raw output of GTSM-ERA5 model runs into daily maximum timeseries of water levels.
"""

import xarray as xr
import sys
import os
import numpy as np
from datetime import datetime

def resampleTS_DailyMax(year):

    print(year)
    year=int(year)

    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_out_wl = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-dailymax','waterlevel') 
    ofile_wl = f'{dir_out_wl}/reanalysis_waterlevel_dailymax_{year}.nc' #output file

    # get list of stations that are included in the previous dataset on CDS
    print('Loading dataset with a list of stations...')
    file_nc_cds = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/reanalysis_waterlevel_hourly_1979_01_v1.nc')
    ds_cds = xr.open_dataset(file_nc_cds); ds_cds.close()
    
    for name in {'waterlevel'}:
    
        # skip if already processed
        if (name=='waterlevel') & (os.path.isfile(ofile_wl)==True):
            print (f'Year {year}, {name} timeseries already processed.')
            continue    

        for mnth in range(1,13):

            # get timeseries
            print(f'Loading timeseries of {name}, month {mnth}...')
            if ((year < 1980) | (year > 2018)):
                file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/{name}/reanalysis_{name}_10min_{year}_{mnth:02d}_v1.nc')
            else:
                file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5-10min-1979-2018/{name}/reanalysis_waterlevel_10min_{year}_{mnth:02d}_v1.nc')
            
            ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
            ds.load()

            #set coordinates
            ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))       
            ds = ds.rename({'station_x_coordinate':'lon','station_y_coordinate':'lat'})

            # remove stations
            print('removing stations...')
            unique_stations = np.setxor1d(ds['stations'], ds_cds['stations'])
            ds = ds.drop_sel(stations=unique_stations)            

            # update attributes
            ds.stations.attrs = {'units': '', 'long_name': 'station_id'}
            ds.waterlevel.attrs = {'long_name': 'Daily maximum of water level above mean sea level','units': 'm','standard_name': 'sea_surface_height_above_mean_sea_level',
 'description': 'Daily maximum of total water level resulting from the combination of barotropic tides and surges and mean sea-level'}
            
            # update general attributes
            ds.attrs={'Conventions':'CF-1.8', 
              'history': 'This is version 1 of the dataset',
              'institution': 'Deltares', 
              'source': 'GTSMv3 forced with ERA5 climate reanalysis',
              'comment':'',
              'references':'DOI: 10.5281/zenodo.10671284',
              'featureType': 'point', 
              'id': 'GTSM-ERA5-E_water_level_daily_maxima',
              'naming_authority': 'https://deltares.nl/en',
              'summary': 'This dataset has been produced with the Global Tide and Surge Model (GTSM) version 3.0. GTSM was forced with wind speed and pressure fields from ERA5 climate reanalysis',
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
                'geospatial_vertical_min': ds[name].load().min().round(3).astype(str).item(), 
                'geospatial_vertical_max': ds[name].load().max().round(3).astype(str).item(),
                'geospatial_vertical_units': 'm', 
                'geospatial_vertical_positive': 'up'}

            # resample dataset to daily max values
            time = ds['time'].values
            time = time[0:-1:144]
            ds = ds.coarsen(time=144).max()
            ds['time'] = time

            if mnth==1:
                ds_all = ds; del ds
            else:
                ds_all = xr.concat([ds_all,ds],dim="time"); del ds
                

        # save daily max timeseries
        ds_all.attrs['time_coverage_start'] = str(ds_all.time.min().dt.strftime('%Y-%m-%d %H:%M:%S').item())   
        ds_all.attrs['time_coverage_end'] = str(ds_all.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())        
        if name=='waterlevel':
            ds_all.attrs['title'] = 'Daily maximum timeseries of total water levels'
            ds_all.to_netcdf(ofile_wl)
    return

if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>0:
        year = int(sys.argv[1])      
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year')
    resampleTS_DailyMax(year)
    
    
    
    
    
    

