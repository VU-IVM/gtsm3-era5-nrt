# Author: Robyn S. Gwee, Sanne Muis, Natalia Aleksandrova
# Contact: robyn.gwee@deltares.nl
# Date created: Mon Feb 22 2021
# Date updated: Thu Sep 14 2023
# Remarks: This script processes the outputs of GTSM-ERA5 model runs into timeseries files of still water levels and surge levels per year and month.

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
        
    # # plot WORLD
    # proj = ccrs.PlateCarree()
    # fig, axes = plt.subplots(ncols=1, nrows=3,figsize=(20,7), subplot_kw={'projection': proj})
    # sc0=axes[0].scatter(x,y,s=10,c=v0,transform=proj,vmin=-0.5,vmax=0.5)
    # fig.colorbar(sc0,ax=axes[0])
    # axes[0].title.set_text('mean (m)')
    # sc1=axes[1].scatter(x,y,s=10,c=v1,transform=proj,vmin=-2,vmax=0)
    # fig.colorbar(sc1,ax=axes[1])
    # axes[1].title.set_text('min (m)')
    # sc2=axes[2].scatter(x,y,s=10,c=v2,transform=proj,vmin=0,vmax=2)
    # fig.colorbar(sc2,ax=axes[2])
    # axes[2].title.set_text('max (m)')
    # for ax in axes.flat:
    #     ax.set_global()
    #     ax.coastlines()
    # fig.savefig(outpath.replace('.nc','.png'))
    
    # plot NL
    # fig, axes = plt.subplots(ncols=1, nrows=3,figsize=(20,7), subplot_kw={'projection': proj})
    # sc0=axes[0].scatter(x,y,s=10,c=v0,transform=proj,vmin=-0.5,vmax=0.5)
    # fig.colorbar(sc0,ax=axes[0])
    # axes[0].title.set_text('mean (m)')
    # sc1=axes[1].scatter(x,y,s=10,c=v1,transform=proj,vmin=-2,vmax=0)
    # fig.colorbar(sc1,ax=axes[1])
    # axes[1].title.set_text('min (m)')
    # sc2 = axes[2].scatter(x,y,s=10,c=v2,transform=proj,vmin=0,vmax=2)
    # fig.colorbar(sc2,ax=axes[2])
    # axes[2].title.set_text('max (m)')
    # for ax in axes.flat:
    #     ax.set_extent([0, 10, 45, 55], ccrs.PlateCarree())
    #     ax.coastlines()
    # fig.savefig(outpath.replace('.nc','_NL.png'))


def exportTWL(dataset, outpath, raw_data):
    print(dataset)
    dataset.waterlevel.attrs = {'long_name': 'sea_surface_height_above_mean_sea_level',
                              'units': 'm',
                              'short_name': 'waterlevel',
                              'description': 'Total water level resulting from the combination of barotropic tides and surges and mean sea-level'}
    dataset.station_y_coordinate.attrs = {'units': 'degrees_north',
                              'short_name': 'latitude',
                              'long_name': 'latitude'}
    dataset.station_x_coordinate.attrs = {'units': 'degrees_east',
                               'short_name': 'longitude',
                                'long_name': 'longitude'}
    dataset.time.attrs = {'axis': 'T',
                            'long_name': 'time',
                            'short_name': 'time'}
    dataset['waterlevel'] = dataset.waterlevel.round(3)
    #dataset = dataset.assign_coords({'stations': dataset.stations})
    encoding={'stations':{'dtype': 'uint16', 'complevel': 3, 'zlib': True},
                'station_y_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                'station_x_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                'waterlevel': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True}}
    dataset=attrib(dataset)
    dataset.attrs['geospatial_vertical_min'] = dataset.waterlevel.min().round(3).astype(str).item()
    dataset.attrs['geospatial_vertical_max'] = dataset.waterlevel.max().round(3).astype(str).item()
    dataset.attrs['source'] = 'GTSMv3 forced with ERA5 climate reanalysis'
    dataset.attrs['id'] = 'GTSMv3_totalwaterlevel'; dataset.attrs['title'] = '10-minute timeseries of total water levels'
    dataset.to_netcdf(outpath, encoding=encoding)
    return dataset


def exportSurge(water, tides, outpath):
    surge_arr = water.waterlevel - tides.tide
    dss = water
    dss = dss.drop(['waterlevel'])
    dss['surge'] = surge_arr.round(3)    
    dss.to_netcdf(outpath.replace('.nc','_tmp.nc'))
    dss = xr.open_dataset(outpath.replace('.nc','_tmp.nc'))
    os.remove(outpath.replace('.nc','_tmp.nc'))
    dss.surge.attrs = {'long_name': 'storm_surge',
                        'units': 'm',
                        'short_name': 'storm_surge',
                        'description': 'Surge signal resulting from subtracting water level and tide variables.'}   
    dss.attrs.update({'id': 'GTSMv3_surge',
                        'title': '10-minute timeseries of surge levels',
                        'source': 'GTSMv3 forced with ERA5 climate reanalysis',
                        'geospatial_vertical_min': dss.surge.min().round(3).astype(str).item(),
                        'geospatial_vertical_max': dss.surge.max().round(3).astype(str).item()
                        })
    
    encoding={'stations':{'dtype': 'uint16', 'complevel': 3, 'zlib': True},
                'station_y_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                'station_x_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
                'surge': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True}}  
    dss.to_netcdf(outpath, encoding=encoding)
    return


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
                  'history': 'This is version 2 of the dataset, that includes the extension of the dataset to the period of 1950-2024. Data for 1980-2018 remains unchanged.', 
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
    raw_data = {'era5': {'fpath': os.path.join(mpath, f'model_input_ERA5_{year}', 'output', 'gtsm_fine_0000_his.nc'),
                        'fname': 'era5', 
                        'opath': 'timeseries-GTSM-ERA5',
                        'opath_stats': 'stats-GTSM-ERA5'},
                }
    print(raw_data)
    inpath = raw_data[scenario]['fpath']
    
    ppath = path_dict['postproc']
    outpath = os.path.join(ppath, raw_data[scenario]['opath'])
    
    wpath = os.path.join(outpath, 'waterlevel')
    spath = os.path.join(outpath, 'surge')  
    
    outpath_stats = os.path.join(ppath, raw_data[scenario]['opath_stats'])
    wpath_stats = os.path.join(outpath_stats, 'waterlevel')
    spath_stats = os.path.join(outpath_stats, 'surge') 
    
    for ipath in [wpath, spath]:
        #if not os.path.exists(ipath):
        os.makedirs(ipath, exist_ok=True)
    del ipath
                
    for ipath in [wpath_stats, spath_stats]:
        #if not os.path.exists(ipath):
        os.makedirs(ipath, exist_ok=True)
    del ipath
            
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
    
    # tidal dataset
    if year>2014:
      texp = 'future'
    else:
      texp = 'historical'
    exp = 'reanalysis'
    
    if (year==1950) & (mnth==1): 
        tpath = os.path.join(outpath, 'tide')
        tfile = os.path.join(tpath,'{}_tide_{}_{}_v2.nc'.format(texp,date.strftime('%Y'), date.strftime('%m')))  
    else:
        tpath = path_dict['tides_CDS']
        tfile = os.path.join(tpath,'{}_tide_{}_{}_v1.nc'.format(texp,date.strftime('%Y'), date.strftime('%m')))  

    # file paths 
    sfile = os.path.join(spath,'{}_surge_10min_{}_{}_v1.nc'.format(exp,date.strftime('%Y'), date.strftime('%m')))
    wfile = os.path.join(wpath,'{}_waterlevel_10min_{}_{}_v1.nc'.format(exp,date.strftime('%Y'),date.strftime('%m')))
    
    if not os.path.exists(wfile):
        print('writing: ', wfile)
        exportTWL(mth, wfile, raw_data)
    resample_and_plot(wfile)
    del mth
    if not os.path.exists(sfile):
        tf = xr.open_dataset(tfile, chunks={'stations': 1000}); tf.close() 
        wf = xr.open_dataset(wfile, chunks={'stations': 1000}); wf.close()
        print('writing: ', sfile)
        exportSurge(wf, tf, sfile)
    resample_and_plot(sfile)
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