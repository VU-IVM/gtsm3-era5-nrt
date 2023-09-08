# -*- coding: utf-8 -*-
"""
Created on Fri  31 2023

@author: n-aleksandrova,veenstra
"""

import xarray as xr
import sys
import os


def resampleTS(year, mnth):

    print(year)
    year=int(year)
    mnth=int(mnth)

    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_out_wl = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly','waterlevel') 
    dir_out_surge = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly','surge')

    # resample TS for water levels
    file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/waterlevel/era5_reanalysis_waterlevel_{year}_{mnth:02d}_v1.nc')
    
    ds = xr.open_dataset(file_nc); ds.close()
    #ds_hr = ds.resample(time='H').mean(dim='time')
    ds_hr = ds.coarsen(time=6).mean()
    del ds
    ds_hr.attrs['title'] = '1-hour timeseries of total water levels'
    ds_hr.attrs['id'] = 'GTSMv3_totalwaterlevel_resampled_1hour'
    print(ds_hr)
    ds_hr.to_netcdf(f'{dir_out_wl}/era5_reanalysis_waterlevel_{year}_{mnth:02d}.nc')
    del ds_hr

    # resample TS for surge
    file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/surge/era5_reanalysis_surge_{year}_{mnth:02d}_v1.nc')
    
    ds = xr.open_dataset(file_nc); ds.close()
    #ds_hr = ds.resample(time='H').mean(dim='time')
    ds_hr = ds.coarsen(time=6).mean()
    del ds
    ds_hr.attrs['title'] = '1-hour timeseries of surge levels'
    ds_hr.attrs['id'] = 'GTSMv3_surge_resampled_1hour'
    ds_hr.to_netcdf(f'{dir_out_surge}/era5_reanalysis_surge_{year}_{mnth:02d}.nc')
    del ds_hr
    return

if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        year = int(sys.argv[1])
        mnth = int(sys.argv[2])        
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year, second - month')
    resampleTS(year, mnth)
    

