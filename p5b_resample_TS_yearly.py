# -*- coding: utf-8 -*-
"""
Created on Fri  31 2023

@author: n-aleksandrova,veenstra
"""

import xarray as xr
import sys
import os


def resampleTS_year(year):

    print(year)
    year=int(year)

    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_out_wl = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly-year-files','waterlevel') 
    dir_out_surge = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly-year-files','surge')

    for mnth in range(1,13):  
        # resample TS for water levels
        file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/waterlevel/era5_reanalysis_waterlevel_{year}_{mnth:02d}_v1.nc')
    
        ds = xr.open_dataset(file_nc); ds.close()
        #ds_hr = ds.resample(time='H').mean(dim='time')
        ds_hr = ds.coarsen(time=6).mean()
        del ds
        if mnth == 1:
            ds_hr_all = ds_hr
        else:
            ds_hr_all = xr.concat([ds_hr_all,ds_hr],dim="time")
        del ds_hr
    
    ds_hr_all.attrs['title'] = '1-hour timeseries of total water levels'
    ds_hr_all.attrs['id'] = 'GTSMv3_totalwaterlevel_resampled_1hour'
    ds_hr_all.attrs['time_coverage_start'] = str(ds_hr_all.time.min().dt.strftime('%Y-%m-%d %H:%M:%S').item())
    ds_hr_all.attrs['time_coverage_end'] = str(ds_hr_all.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())
    print(ds_hr_all)
    ds_hr_all.to_netcdf(f'{dir_out_wl}/era5_reanalysis_waterlevel_{year}.nc')
    del ds_hr_all

    # resample TS for surge
    for mnth in range(1,13):
        file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/surge/era5_reanalysis_surge_{year}_{mnth:02d}_v1.nc')
    
        ds = xr.open_dataset(file_nc); ds.close()
        #ds_hr = ds.resample(time='H').mean(dim='time')
        ds_hr = ds.coarsen(time=6).mean()
        del ds
        if mnth == 1:
            ds_hr_all = ds_hr
        else:
            ds_hr_all = xr.concat([ds_hr_all,ds_hr],dim="time")
        del ds_hr

    ds_hr_all.attrs['title'] = '1-hour timeseries of surge levels'
    ds_hr_all.attrs['id'] = 'GTSMv3_surge_resampled_1hour'
    ds_hr_all.attrs['time_coverage_start'] = str(ds_hr_all.time.min().dt.strftime('%Y-%m-%d %H:%M:%S').item())
    ds_hr_all.attrs['time_coverage_end'] = str(ds_hr_all.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())
    ds_hr_all.to_netcdf(f'{dir_out_surge}/era5_reanalysis_surge_{year}.nc')
    del ds_hr_all
    return

if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>0:
        year = int(sys.argv[1])      
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year')
    resampleTS_year(year)
    

