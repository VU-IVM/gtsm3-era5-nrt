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
  import matplotlib
  # load data
  outpath_stats = outpath.replace('timeseries','stats')
  ds = xr.open_dataset(outpath)
  # compute monthly and annual min, max, mean  
  ds_mean = ds.mean()
  ds_mean.to_netcdf(outpath_stats.replace('.nc','_mean.nc'))
  ds_min = ds.min()
  ds_min.to_netcdf(outpath_stats.replace('.nc','_min.nc'))
  ds_max = ds.max()
  ds_max.to_netcdf(outpath_stats.replace('.nc','_max.nc'))
  # plot
  proj = ccrs.PlateCarree()
  fig, axes = plt.subplots(ncols=0, nrows=3,figsize=(20,7), subplot_kw={'projection': proj}
  ds_mean.plot(ax=axes[0],transform=proj,cmap="magma",robust=True)
  ds_min.plot(ax=axes[1],transform=proj,cmap="magma",robust=True)
  ds_max.plot(ax=axes[2],transform=proj,cmap="magma",robust=True)
  for ax in axes.flat:
    ax.coastlines()
    ax.gridlines()
    gl = ax.gridlines(crs=proj, linewidth=1, color='black', alpha=0.2, linestyle="--")
    gl.ylocator = mticker.FixedLocator(np.arange(-90,90,20))
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 25)) 
    gl.xlabels_top = False
    gl.xlabels_bottom = True
    gl.ylabels_left = True
    gl.ylabels_right = False
  fig.savefig(outpath.replace('nc','.png'))

def exportTWL(dataset, outpath):

    dataset = dataset.drop(['station_id', 'bedlevel', 'timestep', 'station_name'])
    if 'x_velocity' in dataset.variables:
        dataset = dataset.drop(['x_velocity', 'y_velocity'])
        
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
       
    dataset.attrs['source'] = 'GTSMv3 forced with {} dataset'.format(raw_data[scenario]['fname'].split('_')[0])

    dataset.attrs['id'] = 'GTSMv3_totalwaterlevel'; dataset.attrs['title'] = '10-minute timeseries of total water levels'
    dataset.to_netcdf(outpath, encoding=encoding)
    return dataset

def exportSurge(water, tides, outpath):
    surge_arr = water.waterlevel - tides.tide

    dss = water
    dss = dss.drop(['waterlevel'])

    dss['surge'] = surge_arr.round(3)
    
    dss.surge.attrs = {'long_name': 'storm_surge',
                      'units': 'm',
                      'short_name': 'storm_surge',
                      'description': 'Surge signal resulting from subtracting water level and tide variables.'}
    
    dss.attrs.update({'id': 'GTSMv3_surge',
                      'title': '10-minute timeseries of surge levels',
                      'source': 'GTSMv3 forced with {} dataset'.format(raw_data[scenario]['fname'].split('_')[0]),
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

def raw2nc(year, scenario):
  # folder paths
  raw_data = {'era5': {'fpath': '''os.path.join(mpath, 'slr_tide_surge',
                                            'model_input_ERA5_{}'.format(str(year)),
                                            'output', 'gtsm_fine_0000_his.nc')''',
                      'fname': 'era5', 
                      'opath': 'timeseries-GTSM-ERA5',
                      'ts':datetime(2019,1,1,0,0,0),
                      'te':datetime(2022,1,1,0,0,0)},}
  inpath = eval(raw_data[scenario]['fpath'])
    outpath = os.path.join('//projects/0/einf3499/model_runs/', raw_data[scenario]['opath']) 
  tpath = os.path.join(r'//projects/0/einf3499/tides_CDS/')
  wpath = os.path.join(outpath, 'waterlevel')
  spath = os.path.join(outpath, 'surge')  
  
  texp = 'future'
  exp = 'reanalysis'

  for i in [tpath, wpath, spath]:
    if not os.path.exists(i):
      os.makedirs(i)
  ds = xr.open_dataset(inpath); ds.close()
    
  for i in range(1,13):
    date = datetime(year, i, 1)
    
    tfile = os.path.join(tpath, '{}_tide_{}_{}_v1.nc'.format(exp,date.strftime('%Y'), date.strftime('%m')))   
          
    print('Processing {}'.format(date.strftime('%Y-%m')))
    mth = ds.sel(time='{}-{}'.format(date.strftime('%Y'), date.strftime('%m')))

    
    tfile = os.path.join(tpath, 
                         '{}_tide_{}_{}_v1.nc'.format(texp,
                                                      date.strftime('%Y'), 
                                                      date.strftime('%m')))   
    sfile = os.path.join(spath, 
                         '{}_{}_surge_{}_{}_v1.nc'.format(raw_data[scenario]['fname'].split('_')[0],
                                                          exp,
                                                          date.strftime('%Y'), 
                                                          date.strftime('%m')))
    wfile = os.path.join(wpath, 
                         '{}_{}_waterlevel_{}_{}_v1.nc'.format(raw_data[scenario]['fname'].split('_')[0],
                                                               exp,
                                                               date.strftime('%Y'), 
                                                               date.strftime('%m')))

    print('Processing {}'.format(date.strftime('%Y-%m')))
    mth = ds.sel(time='{}-{}'.format(date.strftime('%Y'), date.strftime('%m')))

    if not os.path.exists(wfile):
        exportTWL(mth, wfile)

    if not os.path.exists(sfile):
        tf = xr.open_dataset(tfile); tf.close()
        wf = xr.open_dataset(wfile); wf.close()
        exportSurge(wf, tf, sfile)
        
    else:
        continue
  return

if __name__ == "__main__":
  # import modules
  import os
  # read input arguments
  if len(os.sys.argv)>0:
    year = int(sys.argv[1])
    scenario = sys.argv[2]
    mpath = osys.argv[3]   
  else:
    raise RuntimeError('No arguments were provided\nFirst argument should indicate year.\n Second argument for scenario, third for folder path')
  raw2nc(year, scenario, mpath)
  
  # steps
  # remove spinup
  # compute surge residulas
  # change attributes
  # plot maps
  