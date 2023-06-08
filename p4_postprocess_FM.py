# Author: Robyn S. Gwee, Sanne Muis
# Contact: robyn.gwee@deltares.nl
# Date created: Mon Feb 22 17:25:23 2021
# Remarks: raw2nc 

from datetime import datetime
from dateutil.relativedelta import *
import xarray as xr
import numpy as np
import sys
import glob

def resample_and_plot(outpath):
  # import modules
  import cartopy.crs as ccrs
  import matplotlib.pyplot as plt
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
  varname =  list(ds.keys()) 
  print(ds_min)
  try:
    v0 = ds_mean.surge.values
    v1 = ds_min.surge.values
    v2 = ds_max.surge.values
  except AttributeError:
    v0 = ds_mean.waterlevel.values
    v1 = ds_min.waterlevel.values
    v2 = ds_max.waterlevel.values
  # plot WORLD
  proj = ccrs.PlateCarree()
  fig, axes = plt.subplots(ncols=1, nrows=3,figsize=(20,7), subplot_kw={'projection': proj})
  sc0=axes[0].scatter(x,y,s=10,c=v0,transform=proj,vmin=-1,vmax=0)
  cb0 = fig.colorbar(sc0,ax=axes[0])
  axes[0].title.set_text('mean (m)')
  sc1=axes[1].scatter(x,y,s=10,c=v1,transform=proj,vmin=-0.5,vmax=0.5)
  cb1 = fig.colorbar(sc1,ax=axes[1])
  axes[1].title.set_text('min (m)')
  sc2=axes[2].scatter(x,y,s=10,c=v2,transform=proj,vmin=0,vmax=2)
  cb2 = fig.colorbar(sc2,ax=axes[2])
  axes[2].title.set_text('max (m)')
  for ax in axes.flat:
    ax.set_global()
    ax.coastlines()
  fig.savefig(outpath.replace('.nc','.png'))
  # plot NL
  fig, axes = plt.subplots(ncols=1, nrows=3,figsize=(20,7), subplot_kw={'projection': proj})
  sc0=axes[0].scatter(x,y,s=10,c=v0,transform=proj,vmin=-1,vmax=0)
  cb0 = fig.colorbar(sc0,ax=axes[0])
  axes[0].title.set_text('mean (m)')
  sc1=axes[1].scatter(x,y,s=10,c=v1,transform=proj,vmin=-0.5,vmax=0.5)
  cb1 = fig.colorbar(sc1,ax=axes[1])
  axes[1].title.set_text('min (m)')
  sc2 = axes[2].scatter(x,y,s=10,c=v2,transform=proj,vmin=0,vmax=2)
  cb2 = fig.colorbar(sc2,ax=axes[2])
  axes[2].title.set_text('max (m)')
  for ax in axes.flat:
    ax.set_extent([0, 10, 45, 55], ccrs.PlateCarree())
    ax.coastlines()
  fig.savefig(outpath.replace('.nc','_NL.png'))

def exportTWL(dataset, outpath, raw_data, scenario):
  keys = list(dataset.keys()) # remove all variables except water level
  keys.remove('waterlevel')
  dataset = dataset.drop(keys)
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
  dataset = dataset.assign_coords({'stations': dataset.stations})
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
                'id': 'GTSMv3_tides_totalwaterlevel', 
                'naming_authority': 'https://deltares.nl/en', 
                'Metadata_Conventions': 'Unidata Dataset Discovery v1.0', 
                'title': '10-minute timeseries of total water levels', 
                'summary': 'This dataset has been produced with the Global Tide and \
                            Surge Model (GTSM) version 3.0. GTSM was forced with \
                            wind speed and pressure fields from ERA5 climate reanalysis', 
                'date_created': str(datetime.utcnow()) + ' UTC', 
                'project': 'GTSMip and C3S_435_Lot8 Deltares', 
                'acknowledgment': 'The development of this dataset was financed with Deltares Strategic Research Program. Additional funding was received by Contract C3S_435_Lot8 Deltares', 
                'contact': 'kun.yan@deltares.nl',
                'license': 'Copernicus Products License', 
                'history': '', 
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

def raw2nc(scenario,year,mnth, mpath):
  raw_data = {'era5': {'fpath': '''os.path.join(mpath,'model_input_ERA5_{}'.format(str(year)),
                                            'output', 'gtsm_fine_0000_his.nc')''',
                      'fname': 'era5', 
                      'opath': 'timeseries-GTSM-ERA5',
                      'ts':datetime(2019,1,1,0,0,0),
                      'te':datetime(2022,1,1,0,0,0)},}
  print(raw_data)
  inpath = eval(raw_data[scenario]['fpath'])
  outpath = os.path.join(r'/projects/0/einf3499/model_runs/', raw_data[scenario]['opath']) 
  tpath = os.path.join(r'/projects/0/einf3499/tides_CDS/')
  wpath = os.path.join(outpath, 'waterlevel')
  spath = os.path.join(outpath, 'surge')  
  
  texp = 'future'
  exp = 'reanalysis'
  
  for i in [tpath, wpath, spath]:
    if not os.path.exists(i):
      os.makedirs(i)
  ds = xr.open_dataset(inpath); ds.close()
  print('ds loaded')   
  
  date = datetime(year, mnth , 1)      
  print('Processing {}'.format(date.strftime('%Y-%m')))
  mth = ds.sel(time='{}-{}'.format(date.strftime('%Y'), date.strftime('%m')))
  
  tfile = os.path.join(tpath,'{}_tide_{}_{}_v1.nc'.format(texp,date.strftime('%Y'), date.strftime('%m')))   
  sfile = os.path.join(spath,'{}_{}_surge_{}_{}_v1.nc'.format(raw_data[scenario]['fname'].split('_')[0],exp,date.strftime('%Y'), date.strftime('%m')))
  wfile = os.path.join(wpath,'{}_{}_waterlevel_{}_{}_v1.nc'.format(raw_data[scenario]['fname'].split('_')[0],exp,date.strftime('%Y'),date.strftime('%m')))
  
  if not os.path.exists(wfile):
    print('writing: ', wfile)
    exportTWL(mth, wfile, raw_data, scenario)
  resample_and_plot(wfile)
  if not os.path.exists(sfile):
    tf = xr.open_dataset(tfile, chunks={'stations': 1000}); tf.close()
    wf = xr.open_dataset(wfile, chunks={'stations': 1000}); wf.close()
    print('writing: ', sfile)
    exportSurge(wf, tf, sfile)
  resample_and_plot(sfile)
  return
      
if __name__ == "__main__":    
  # import modulesyear=2019
  import os
  # read input arguments
  print(len(os.sys.argv))
  if len(os.sys.argv)>1:
    year = int(sys.argv[1])
    mnth = int(sys.argv[2])
    scenario = sys.argv[3]
    mpath = sys.argv[4]   
  else:
    raise RuntimeError('No arguments were provided\nFirst argument should indicate year.\n Second argument for scenario, third for folder path')
  
  print(year, scenario, mnth, mpath)
  raw2nc(scenario,year, mnth, mpath)