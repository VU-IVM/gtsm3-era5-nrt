# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 2023

@author: n-aleksandrova,veenstra
"""

import xarray as xr
import sys
import os
year = int(sys.argv[1])
print(year)

sys.path.append("..")
from path_dict import path_dict
dir_postproc = path_dict['postproc']
dir_output = os.path.join(dir_postproc,'NL_stations','timeseries') 

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

station_list = ['NWS_NO_TS_MO_Vlissingen','NWS_NO_TS_MO_HoekVanHolland','NWS_NO_TS_MO_Ijmuiden','NWS_NO_TS_MO_DenHelder','NWS_NO_TS_MO_Harlingen','NWS_NO_TS_MO_Delfzijl']

# get timeseries for water levels
for mnth in range(1,13):
    file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/waterlevel/era5_reanalysis_waterlevel_{year}_{mnth:02d}_v1.nc')
    ds = xr.open_dataset(file_nc,chunks={'stations':1000}); ds.close()
    
    #set station_name variable/coordinate for stations dimensions as its index
    ds['station_name'] = ds['station_name'].load().str.decode('utf-8',errors='ignore').str.strip()
    ds = ds.set_index({'stations':'station_name'})

    #drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
    duplicated_keepfirst = ds['stations'].to_series().duplicated(keep='first')
    if duplicated_keepfirst.sum()>0:
        print(f'dropping {duplicated_keepfirst.sum()} duplicate "station_name" labels to avoid InvalidIndexError')
        ds = ds[{'stations':~duplicated_keepfirst}]
    
    # select specific stations        
    data_sel = ds.sel(stations=station_list)
    
    if mnth == 1:
        data_sel_all = data_sel
    else:
        data_sel_all = xr.concat([data_sel_all,data_sel],dim="time")
print(data_sel_all)

# resample to 1 hour
data_hourly = data_sel_all.resample(time='H').mean(dim='time')
data_hourly.attrs['title'] = '1-hour timeseries of total water levels'
data_hourly.attrs['time_coverage_start'] = str(data_hourly.time.min().dt.strftime('%Y-%m-%d %H:%M:%S').item())
data_hourly.attrs['time_coverage_end'] = str(data_hourly.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())
data_hourly.attrs['id'] = 'GTSMv3_totalwaterlevel_resampled_1hour'
print(data_hourly)

data_hourly.to_netcdf(f'{dir_output}/era5_reanalysis_waterlevel_{year}_NLstations.nc')

# get timeseries for surge
for mnth in range(1,13):
    file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/surge/era5_reanalysis_surge_{year}_{mnth:02d}_v1.nc')
    ds = xr.open_dataset(file_nc,chunks={'stations':1000}); ds.close()
    
    #set station_name variable/coordinate for stations dimensions as its index
    ds['station_name'] = ds['station_name'].load().str.decode('utf-8',errors='ignore').str.strip()
    ds = ds.set_index({'stations':'station_name'})

    #drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
    duplicated_keepfirst = ds['stations'].to_series().duplicated(keep='first')
    if duplicated_keepfirst.sum()>0:
        print(f'dropping {duplicated_keepfirst.sum()} duplicate "station_name" labels to avoid InvalidIndexError')
        ds = ds[{'stations':~duplicated_keepfirst}]
    
    # select specific stations        
    data_sel = ds.sel(stations=station_list)
    
    if mnth == 1:
        data_sel_all = data_sel
    else:
        data_sel_all = xr.concat([data_sel_all,data_sel],dim="time")
print(data_sel_all)

# resample to 1 hour
data_hourly = data_sel_all.resample(time='H').mean(dim='time')
data_hourly.attrs['title'] = '1-hour timeseries of surge levels'
data_hourly.attrs['time_coverage_start'] = str(data_hourly.time.min().dt.strftime('%Y-%m-%d %H:%M:%S').item())
data_hourly.attrs['time_coverage_end'] = str(data_hourly.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())
data_hourly.attrs['id'] = 'GTSMv3_surge_resampled_1hour'
print(data_hourly)

data_hourly.to_netcdf(f'{dir_output}/era5_reanalysis_surge_{year}_NLstations.nc')



    

