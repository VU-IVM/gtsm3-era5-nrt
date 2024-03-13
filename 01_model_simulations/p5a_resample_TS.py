# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 2023
@author: n-aleksandrova,veenstra
"""

import xarray as xr
import os
import numpy as np
import sys
sys.path.append("..")
from path_dict import path_dict

def resampleTS(year, mnth):

    print(year)
    year=int(year)
    mnth=int(mnth)

    dir_postproc = path_dict['postproc']
    dir_out_wl = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly','waterlevel') 
    dir_out_surge = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly','surge')
    dir_out_wl_10min = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-10min','waterlevel') 
    dir_out_surge_10min = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-10min','surge')

    ofile_wl_1hr = f'{dir_out_wl}/reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc' #output file
    ofile_surge_1hr = f'{dir_out_surge}/reanalysis_surge_hourly_{year}_{mnth:02d}_v1.nc' #output file
    
    ofile_wl_10min = f'{dir_out_wl_10min}/reanalysis_waterlevel_10min_{year}_{mnth:02d}_v1.nc' #output file
    ofile_surge_10min = f'{dir_out_surge_10min}/reanalysis_surge_10min_{year}_{mnth:02d}_v1.nc' #output file
    
    # get list of stations that are included in the previous dataset on CDS
    print('Loading dataset with a list of stations...')
    file_nc_cds = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/reanalysis_waterlevel_hourly_1979_01_v1.nc')
    ds_cds = xr.open_dataset(file_nc_cds); ds_cds.close()
    
    for name in {'waterlevel','surge'}:
    
        # skip if already processed
        if (name=='waterlevel') & (os.path.isfile(ofile_wl_1hr)==True) & (os.path.isfile(ofile_wl_10min)==True):
            print (f'Year {year}, month {mnth}, {name} timeseries already processed.')
            continue

        if (name=='surge') & (os.path.isfile(ofile_surge_1hr)==True) & (os.path.isfile(ofile_surge_10min)==True):
            continue         

        # get timeseries
        print(f'Loading timeseries of {name}...')
        file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/{name}/era5_reanalysis_{name}_10min_{year}_{mnth:02d}_v1.nc')
        ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
        ds.load()

        # remove unnecessary coordinates
        try:
            ds = ds.drop(['station_name'])
        except ValueError:
            print('Key station_name already removed')

        #set coordinates
        ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))       

        # remove stations
        print('removing stations...')
        unique_stations = np.setxor1d(ds['stations'], ds_cds['stations'])
        ds = ds.drop_sel(stations=unique_stations)            
    
        # update general attributes
        ds.attrs['contact'] = '', #'Please contact Copernicus User Support on the Copernicus Climate Change Service website ( https://climate.copernicus.eu/).'
        ds.attrs['geospatial_vertical_min'] = ds[name].load().min().round(3).astype(str).item() 
        ds.attrs['geospatial_vertical_max'] = ds[name].load().max().round(3).astype(str).item() 
        ds.attrs['geospatial_lat_min'] = ds.station_y_coordinate.min().round(3).astype(str).item()
        ds.attrs['geospatial_lat_max'] = ds.station_y_coordinate.max().round(3).astype(str).item() 
        ds.attrs['geospatial_lon_min'] = ds.station_x_coordinate.min().round(3).astype(str).item()
        ds.attrs['geospatial_lon_max'] = ds.station_x_coordinate.max().round(3).astype(str).item()
        ds.attrs['project'] = 'Deltares Strategic Research Program' 
        ds.attrs['acknowledgment'] = 'The development of this dataset was financed with Deltares Strategic Research Program.'

        # save 10-minute timeseries
        print(ds)
        if name=='waterlevel':
            ds.to_netcdf(ofile_wl_10min)
        elif name=='surge':
            ds.to_netcdf(ofile_surge_10min)

        # resample dataset to hourly values
        time = ds['time'].values
        time = time[0:-1:6]
        ds = ds.coarsen(time=6).mean()
        ds['time'] = time

        # save hourly timeseries
        ds.attrs['time_coverage_end'] = str(ds.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())        
        if name=='waterlevel':
            ds.attrs['title'] = 'Hourly timeseries of total water levels'
            ds.to_netcdf(ofile_wl_1hr)
        elif name=='surge':
            ds.attrs['title'] = 'Hourly timeseries of surge levels'
            ds.to_netcdf(ofile_surge_1hr)
    return

if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        year = int(sys.argv[1])
        mnth = int(sys.argv[2])        
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year, second - month')
    resampleTS(year, mnth)
    
    
    
    
    
    

