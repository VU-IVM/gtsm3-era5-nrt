# Author: Robyn S. Gwee, Sanne Muis, Natalia Aleksandrova
# Date created: Mon Feb 22 2021
# Date updated: Thu May 6 2025
# Remarks: This script processes the outputs of GTSM-ERA5 tide-only model runs into timeseries files 

import os
from datetime import datetime
import xarray as xr
import sys
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import sys
import numpy as np
sys.path.append("..")
from path_dict import path_dict


def resample_and_plot(outpath):
    # load data
    outpath_stats = outpath.replace('timeseries','stats')
    ds = xr.open_dataset(outpath)
    print(ds.keys())
    # compute monthly min, max, mean  
    ds_mean = ds.mean(dim='time')
    ds_mean.to_netcdf(outpath_stats.replace('.nc','_monthly_mean.nc'))
    #ds_mean = ds_mean.sel(station_y_coordinate=slice(45,55), station_x_coordinate=slice(0,10))
    #ds_mean.to_netcdf(outpath_stats.replace('.nc','_monthly_mean_NL.nc'))
    ds_min = ds.min(dim='time')
    ds_min.to_netcdf(outpath_stats.replace('.nc','_monthly_min.nc'))
    ds_max = ds.max(dim='time')
    ds_max.to_netcdf(outpath_stats.replace('.nc','_monthly_max.nc')) 
    # coor
    x = ds['station_x_coordinate']
    y = ds['station_y_coordinate']
    print(ds_min)
    try:
        v0 = ds_mean.surge.values
        v1 = ds_min.surge.values
        v2 = ds_max.surge.values
    except AttributeError:
        v0 = ds_mean.waterlevel.values
        v1 = ds_min.waterlevel.values
        v2 = ds_max.waterlevel.values

def exportTide(dataset, outpath, raw_data):
    tide_arr = dataset.waterlevel
    dataset = dataset.drop(['waterlevel'])
    dataset['tide'] = tide_arr.round(3)    
    dataset.to_netcdf(outpath.replace('.nc','_tmp.nc'))
    del dataset
    dataset = xr.open_dataset(outpath.replace('.nc','_tmp.nc'))
    os.remove(outpath.replace('.nc','_tmp.nc'))
    
    print(dataset)
    
    dataset.tide.attrs = {'long_name': 'sea_surface_height_amplitude_due_to_tide',
                              'units': 'm',
                              'short_name': 'tide',
                              'description': 'Barotropic tidal signal containing astronomic tide, self-attraction and loading, radiational tides and mean sea-level'}
    dataset.station_y_coordinate.attrs = {'units': 'degrees_north',
                              'short_name': 'latitude',
                              'long_name': 'latitude'}
    dataset.station_x_coordinate.attrs = {'units': 'degrees_east',
                               'short_name': 'longitude',
                                'long_name': 'longitude'}
    dataset.time.attrs = {'axis': 'T',
                            'long_name': 'time',
                            'short_name': 'time'}
    #dataset = dataset.assign_coords({'stations': dataset.stations})
    encoding={'stations':{'dtype': 'uint16', 'complevel': 3, 'zlib': True},
                'station_y_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                'station_x_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                'tide': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True}}
    dataset=attrib(dataset)
    dataset.attrs['geospatial_vertical_min'] = dataset.tide.min().round(3).astype(str).item()
    dataset.attrs['geospatial_vertical_max'] = dataset.tide.max().round(3).astype(str).item()
    dataset.attrs['source'] = 'GTSMv3'
    dataset.attrs['id'] = 'GTSMv3_tides'; dataset.attrs['title'] = '10-minute timeseries of tide levels'
    dataset.to_netcdf(outpath, encoding=encoding)
    return dataset

def attrib(dataset):
    # add global attribute metadata
    if 'tide' in list(dataset.keys()):
        gvm = 'tide'
    elif 'waterlevel' in list(dataset.keys()):
        gvm = 'waterlevel'
    elif 'surge' in list(dataset.keys()):
        gvm = 'surge'
    experiment = 'reanalysis'
    
    dataset.attrs={'Conventions':'CF-1.6', 
                  'featureType': 'timeSeries', 
                  'id': 'GTSMv3_totalwaterlevels', 
                  'naming_authority': 'https://deltares.nl/en', 
                  'Metadata_Conventions': 'Unidata Dataset Discovery v1.0', 
                  'title': '10-minute timeseries of total water levels', 
                  'summary': 'This dataset has been produced with the Global Tide and Surge Model (GTSM) version 3.0. GTSM was forced with wind speed and pressure fields from ERA5 climate reanalysis', 
                  'date_created': str(datetime.utcnow()) + ' UTC', 
                  'date_modified': '', 
                  'project': 'Deltares Strategic Research Program', 
                  'acknowledgment': 'The development of this dataset was financed with Deltares Strategic Research Program.', 
                  'contact': '',
                  'license': 'Copernicus Products License', 
                  'history': 'This is version 1 of the dataset', 
                  'institution': 'Deltares', 
                  'sea_name': 'global', 
                  'source': 'GTSMv3 forced with ERA5 climate reanalysis',
                  'keywords': 'sea-level rise; climate change; water level; climate; tides; hydrography; global tide and surge model;', 
                  'keywords_vocabulary': 'http://www.eionet.europa.eu/gemet', 
                  'standard_name_vocabulary': '', 
                  'geospatial_lat_min': dataset.station_y_coordinate.min().round(3).astype(str).item(), 
                  'geospatial_lat_max': dataset.station_y_coordinate.max().round(3).astype(str).item(), 
                  'geospatial_lon_min': dataset.station_x_coordinate.min().round(3).astype(str).item(), 
                  'geospatial_lon_max': dataset.station_x_coordinate.max().round(3).astype(str).item(), 
                  'geospatial_lat_units': 'degrees_north',
                  'geospatial_lat_resolution': 'point',
                  'geospatial_lon_units': 'degrees_east', 
                  'geospatial_lon_resolution': 'point',
                  'geospatial_vertical_min': dataset[gvm].min().round(3).astype(str).item(), 
                  'geospatial_vertical_max': dataset[gvm].max().round(3).astype(str).item(),
                  'geospatial_vertical_units': 'm', 
                  'geospatial_vertical_positive': 'up',
                  'time_coverage_start': str(dataset.time.min().dt.strftime('%Y-%m-%d %H:%M:%S').item()), 
                  'time_coverage_end': str(dataset.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item()),
                  'experiment': experiment}
    return dataset


def raw2nc(year, mnth, scenario):
    mpath = path_dict['modelruns']
    raw_data = {'era5': {'fpath': os.path.join(mpath, f'model_input_TIDE_{year}', 'output', 'gtsm_fine_0000_his.nc'),
                        'fname': 'tide', 
                        'opath': 'timeseries-GTSM-ERA5',
                        'opath_stats': 'stats-GTSM-ERA5'},
                }
    print(raw_data)
    inpath = raw_data[scenario]['fpath']
    
    ppath = path_dict['postproc']
    outpath = os.path.join(ppath, raw_data[scenario]['opath']) 
    tpath = os.path.join(outpath, 'tide')  
    
    outpath_stats = os.path.join(ppath, raw_data[scenario]['opath_stats'])
    tpath_stats = os.path.join(outpath_stats, 'tide') 
    
    if year>2014:
      texp = 'future'
    else:
      texp = 'historical'
    exp = 'reanalysis'
    
    os.makedirs(tpath, exist_ok=True)
    os.makedirs(tpath_stats, exist_ok=True)
            
    ds = xr.open_dataset(inpath);
    print('ds loaded')   
    
    date = datetime(year, mnth , 1)
    print('Processing {}'.format(date.strftime('%Y-%m')))
    mth = ds.sel(time='{}-{}'.format(date.strftime('%Y'), date.strftime('%m')),drop=True)
    mth.load()
    ds.close()

    keys = list(mth.keys()) # remove all variables except water level
    keys.remove('waterlevel')
    mth = mth.drop(keys)
    mth = mth.assign_coords({'stations': mth.stations})

    # subset stations to match stations used in the previous GTSM-ERA5 dataset - filtered stations to remove erroneous points (e.g. too much inland etc.)
    print('Loading dataset with a list of stations...')
    dir_postproc = path_dict['postproc']
    file_nc_cds = os.path.join(dir_postproc,f'timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/reanalysis_waterlevel_hourly_1979_01_v1.nc')
    ds_cds = xr.open_dataset(file_nc_cds); 

    # remove stations
    print('removing stations...')
    unique_stations = np.setxor1d(mth['stations'], ds_cds['stations'])
    mth = mth.drop_sel(stations=unique_stations)   
    ds_cds.close()

    # file paths
    tfile = os.path.join(tpath,'{}_tide_{}_{}_v2.nc'.format(texp,date.strftime('%Y'), date.strftime('%m')))   
    
    if not os.path.exists(tfile):
        print('writing: ', tfile)
        exportTide(mth, tfile, raw_data)
    else:
        print('File already exists: ', tfile)
    #resample_and_plot(tfile)
    del mth
    return


if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        year = int(sys.argv[1])
        mnth = int(sys.argv[2])
        scenario = sys.argv[3]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year.\n Third argument for scenario')
    
    print(year, mnth, scenario)
    raw2nc(year, mnth, scenario)