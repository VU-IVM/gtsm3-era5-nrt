# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 2023
@author: n-aleksandrova,veenstra
"""

import xarray as xr
import sys
import os
import numpy as np

def resampleTS_DailyMax(year):

    print(year)
    year=int(year)

    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_out_wl = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-dailymax','waterlevel') 
    ofile_wl = f'{dir_out_wl}/reanalysis_waterlevel_hourly_{year}.nc' #output file

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
            file_nc = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5/{name}/era5_reanalysis_{name}_10min_{year}_{mnth:02d}_v1.nc')
            ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
            ds.load()

#            # remove unnecessary coordinates
#            try:
#                ds = ds.drop(['station_name'])
#            except ValueError:
#                print('Key station_name already removed')

            #set coordinates
            ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))       

            # remove stations
            print('removing stations...')
            unique_stations = np.setxor1d(ds['stations'], ds_cds['stations'])
            ds = ds.drop_sel(stations=unique_stations)            
    
            # update general attributes
            ds.attrs['contact'] = 'natalia.aleksandrova@deltares.nl', #'Please contact Copernicus User Support on the Copernicus Climate Change Service website ( https://climate.copernicus.eu/).'
            ds.attrs['geospatial_vertical_min'] = ds[name].load().min().round(3).astype(str).item() 
            ds.attrs['geospatial_vertical_max'] = ds[name].load().max().round(3).astype(str).item() 
            ds.attrs['geospatial_lat_min'] = ds.station_y_coordinate.min().round(3).astype(str).item()
            ds.attrs['geospatial_lat_max'] = ds.station_y_coordinate.max().round(3).astype(str).item() 
            ds.attrs['geospatial_lon_min'] = ds.station_x_coordinate.min().round(3).astype(str).item()
            ds.attrs['geospatial_lon_max'] = ds.station_x_coordinate.max().round(3).astype(str).item()
            ds.attrs['project'] = 'Deltares Strategic Research Program' 
            ds.attrs['acknowledgment'] = 'The development of this dataset was financed with Deltares Strategic Research Program.'

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
    resampleTS_DailyMax(year, mnth)
    
    
    
    
    
    

