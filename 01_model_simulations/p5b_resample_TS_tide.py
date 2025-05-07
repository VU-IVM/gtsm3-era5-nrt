# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 2023
@author: n-aleksandrova,veenstra
"""

import xarray as xr
import os
import numpy as np
import sys
import pandas as pd
from datetime import datetime, timezone, timedelta
sys.path.append("..")
from path_dict import path_dict

def resampleTS(year, mnth):

    print(year)
    year=int(year)
    mnth=int(mnth)
    
    # settings for which resampling to apply (by default all)
    save_10min = 1

    dir_postproc = path_dict['postproc']
    dir_out_tide_10min = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-10min','tide') 

    os.makedirs(dir_out_tide_10min, exist_ok=True)

    ofile_tide_10min = f'{dir_out_tide_10min}/historical_tide_{year}_{mnth:02d}_v2.nc' #output file
        
    for name in {'tide'}:

        # get timeseries
        print(f'Loading timeseries of {name}...')
        file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/{name}/historical_{name}_{year}_{mnth:02d}_v2.nc')
        ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
        ds = ds.drop_vars(['station_name'])
        ds.load()

        #set coordinates
        ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))                
    
        # update attributes
        ds.attrs['contact'] = 'Please contact Copernicus User Support on the Copernicus Climate Change Service website (https://climate.copernicus.eu/).'
        ds.attrs['geospatial_vertical_min'] = ds[name].load().min().round(3).astype(str).item() 
        ds.attrs['geospatial_vertical_max'] = ds[name].load().max().round(3).astype(str).item() 
        ds.attrs['project'] = 'Deltares Strategic Research Program' 
        ds.attrs['history'] = 'This is version 2 of the dataset, that includes the extension of the dataset to the period of 1950-2024. Data for 1980-2018 remains unchanged.', 
        ds.attrs['acknowledgment'] = 'The development of this dataset was financed with Deltares Strategic Research Program. Additional funding was received by Contract C3S2_422 Deltares'
        ds.attrs['date_modified'] = str(datetime.now(timezone.utc).replace(tzinfo=None)) + ' UTC'
        ds.attrs['summary'] = 'This dataset has been produced with the Global Tide and Surge Model (GTSM) version 3.0. GTSM was forced with wind speed and pressure fields from ERA5 climate reanalysis.'
        
        ds['time'] = ds.time.assign_attrs({'axis': 'T', 'long_name': 'time', 'short_name': 'time'})
        ds['station_x_coordinate'].attrs['crs'] = 'EPSG:4326'
        ds['station_y_coordinate'].attrs['crs'] = 'EPSG:4326'
        
        encoding={'stations':{'dtype': 'uint16', 'complevel': 3, 'zlib': True},
                'station_y_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                'station_x_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                name: {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True}}
        

        # save 10-minute timeseries
        print(ds)
        if save_10min:
            if (name=='tide') & (os.path.isfile(ofile_tide_10min)==False):
                ds.to_netcdf(ofile_tide_10min, encoding=encoding)
            else:
                print (f'Year {year}, month {mnth}, {name}, 10min timeseries already processed,skipping...')
        
        # resample dataset to hourly values
        if save_hourly: 
           
            # resample to hourly
            ds_offset = ds.copy()
            ds_offset["time"] = ds_offset["time"] + pd.Timedelta(minutes=30)
            ds_hr = ds_offset.resample(time='h').mean()
            
            # adjust attributes
            ds_hr['time'] = ds_hr.time.assign_attrs({'axis': 'T', 'long_name': 'time', 'short_name': 'time'})
            ds_hr.attrs['geospatial_vertical_min'] = ds_hr[name].load().min().round(3).astype(str).item() 
            ds_hr.attrs['geospatial_vertical_max'] = ds_hr[name].load().max().round(3).astype(str).item() 
            ds_hr.attrs['time_coverage_end'] = str(ds_hr.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())        
                
            # save hourly timeseries
            if (name=='waterlevel') & (os.path.isfile(ofile_wl_1hr)==False):
                ds_hr.attrs['title'] = 'Hourly timeseries of total water levels'
                ds_hr.to_netcdf(ofile_wl_1hr, encoding=encoding)
            elif (name=='surge') & (os.path.isfile(ofile_surge_1hr)==False):
                ds_hr.attrs['title'] = 'Hourly timeseries of surge levels'
                ds_hr.to_netcdf(ofile_surge_1hr, encoding=encoding)
            ds_hr.close(); del ds_hr

        # resample dataset to daily max values
        if save_dailymax:
            time = ds['time'].values
            time = time[0:-1:144]
            ds_dmax = ds.coarsen(time=144).max()
            ds_dmax['time'] = time
            ds_dmax['time'] = ds_dmax.time.assign_attrs({'axis': 'T', 'long_name': 'time', 'short_name': 'time'})

            ds_dmax.attrs['geospatial_vertical_min'] = ds_dmax[name].load().min().round(3).astype(str).item() 
            ds_dmax.attrs['geospatial_vertical_max'] = ds_dmax[name].load().max().round(3).astype(str).item() 
            ds_dmax.attrs['time_coverage_start'] = str(ds_dmax.time.min().dt.strftime('%Y-%m-%d %H:%M:%S').item())  
            ds_dmax.attrs['time_coverage_end'] = str(ds_dmax.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())        

            # save hourly timeseries
            if (name=='waterlevel') & (os.path.isfile(ofile_wl_dailymax)==False):
                ds_dmax.attrs['title'] = 'Daily timeseries of total water levels'
                ds_dmax.to_netcdf(ofile_wl_dailymax)
            elif (name=='surge') & (os.path.isfile(ofile_surge_dailymax)==False):
                ds_dmax.attrs['title'] = 'Daily timeseries of surge levels'
                ds_dmax.to_netcdf(ofile_surge_dailymax)
            ds_dmax.close(); del ds_dmax

        ds.close()
        del ds
    return

if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        year = int(sys.argv[1])
        mnth = int(sys.argv[2])        
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year, second - month')
    resampleTS(year, mnth)
    
    
    
    
    
    

