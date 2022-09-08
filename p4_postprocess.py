# Author: Robyn S. Gwee
# Contact: robyn.gwee@deltares.nl
# Date created: Mon Feb 22 17:25:23 2021
# Remarks: raw2nc 
import os
from datetime import datetime
from dateutil.relativedelta import *
import xarray as xr
import numpy as np
import sys
import glob

mpath = os.path.join(r'/p/11205312-006-global_projection_sea/01_RUNS')
#mpath = os.path.join(r'P:\11205312-006-global_projection_sea\01_RUNS')
#tides = os.path.join(r"/p/11205312-006-global_projection_sea/01_RUNS/00_TIDES/tide")
#year = 1960
#scenario = 'gfdl'
year = int(sys.argv[1])
scenario = sys.argv[2]
raw_data = {'tides': {'fpath': '''os.path.join(mpath, '00_TIDES', 'tide', 'model_input_EC-Earth3P-HR_{}'.format(str(year)),
                          'output', 'gtsm_fine_0000_his.nc')''', 
                      'fname': 'tides',
                      'opath': '00_TIDES',
                       'ts':datetime(1950,1,1,0,0,0),
                       'te':datetime(2050,1,1,0,0,0)},
             'hist': {'fpath': '''os.path.join(mpath, '01_HIST', 'slr_tide_surge',
                                            'model_input_EC-Earth3P-HR_{}'.format(str(year)),
                                            'output', 'gtsm_fine_0000_his.nc')''', 
                      'fname': 'EC-Earth3P-HR_hist-1950',
                      'opath': '01_HIST',
                      'ts':datetime(1950,1,1,0,0,0),
                      'te':datetime(2014,12,31,18,10,0)},
             'rcp85': {'fpath': '''os.path.join(mpath, '02_RCP85', 'slr_tide_surge',
                                             'model_input_EC-Earth3P-HR_{}'.format(str(year)),
                                             'output', 'gtsm_fine_0000_his.nc')''',
                       'fname': 'EC-Earth3P-HR_highres-future',
                       'opath': '02_RCP85',
                       'ts':datetime(2015,1,1,0,0,0),
                       'te':datetime(2050,12,31,18,10,0)},
             'era5': {'fpath': '''os.path.join(mpath, '03_ERA5', 'slr_tide_surge',
                                            'model_input_ERA5_{}'.format(str(year)),
                                            'output', 'gtsm_fine_0000_his.nc')''',
                      'fname': 'era5', 
                      'opath': '03_ERA5',
                      'ts':datetime(1979,1,1,0,0,0),
                      'te':datetime(2019,1,1,0,0,0)},
             'sst_hist': {'fpath': '''os.path.join(mpath, '04_SST_HIST', 'slr_tide_surge',
                                                'model_input_SST-HIST_{}'.format(str(year)),
                                                'output', 'gtsm_fine_0000_his.nc')''', 
                          'fname': 'HadGEM3-GC31-HM_highresSST-future',
                          'opath': '04_SST_HIST',
                          'ts': datetime(1951,1,1,0,0,0),
                          'te': datetime(2015,1,1,0,0,0)},
             'sst_rcp85': {'fpath': '''os.path.join(mpath, '05_SST_RCP85', 'slr_tide_surge',
                                                 'model_input_SST-FUT_{}'.format(str(year)),
                                                 'output', 'gtsm_fine_0000_his.nc')''',
                           'fname': 'HadGEM3-GC31-HM_highresSST-present',
                           'opath': '05_SST_RCP85',
                           'ts':datetime(2016,1,1,0,0,0),
                           'te':datetime(2051,1,1,0,0,0)},
             'gfdl': {'fpath':'''glob.glob(os.path.join(mpath, '10_GFDL_HIST', 
                                                      'model_input_GFDL_{}*'.format(str(year)), 
                                                      'output', 'gtsm_fine_0000_his.nc'))[0]''', 
                      'fname': 'GFDL-CM4C192_highresSST-future',
                      'opath': '10_GFDL_HIST',
                      'ts': '1950-01', 
                      'te': '2015-01'},
             'gfdl_rcp85': {'fpath': """glob.glob(os.path.join(mpath, '11_GFDL_RCP85',
                                                            'model_input_GFDL_{}*'.format(str(year)), 
                                                            'output', 'gtsm_fine_0000_his.nc'))[0]""", 
                            'fname': 'GFDL-CM4C192_highresSST-scenario',
                            'opath': '11_GFDL_RCP85'}}

def exportTides(dataset, outpath):

    dataset = dataset.rename({'waterlevel': 'tide'})
    if (year == 1986) or (year == 1985):
        dataset = dataset.drop(['timestep', 'bedlevel'])
    else:
        dataset = dataset.drop(['station_id', 'bedlevel', 'timestep', 'station_name'])
    dataset.tide.attrs = {'long_name': 'sea_surface_height_amplitude_due_to_tide',
                          'units': 'm',
                          'short_name': 'tide',
                          'CDI_grid_type': 'unstructured',
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
    dataset['tide'] = dataset.tide.round(3)
    
    dataset = dataset.assign_coords({'stations': dataset.stations})

    encoding={'stations':{'dtype': 'uint16', 'complevel': 3, 'zlib': True},
              'station_y_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
              'station_x_coordinate': {'dtype': 'int32', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True},
              'tide': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -999, 'complevel': 3, 'zlib': True}}
    
    dataset=attrib(dataset)

    dataset.attrs['id'] = 'GTSMv3_tides'; dataset.attrs['title'] = '10-minute timeseries of tides'
    
    dataset.to_netcdf(outpath, encoding=encoding)
    return

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

def raw2nc(year, scenario):
    inpath = eval(raw_data[scenario]['fpath'])
    outpath = os.path.join('/p/11205028-c3s_435/01_data/01_Timeseries', raw_data[scenario]['opath'],
                           'raw2nc')
    
    tpath = os.path.join(r'/p/11205028-c3s_435/01_data/01_Timeseries/00_TIDES/raw2nc/tide')
    wpath = os.path.join(outpath, 'waterlevel')
    spath = os.path.join(outpath, 'surge')
    if scenario == 'tides':
        if not os.path.exists(tpath):
            os.makedirs(tpath)
    else:
        for i in [tpath, wpath, spath]:
            if not os.path.exists(i):
                os.makedirs(i)
    
    ds = xr.open_dataset(inpath); ds.close()
    
    for i in range(1,13):
        date = datetime(year, i, 1)
        
        if scenario == 'tides':
            if year < 2015:
                exp = 'historical'
            else:
                exp = 'future'
            
            tfile = os.path.join(tpath, 
                             '{}_tide_{}_{}_v1.nc'.format(exp,
                                                          date.strftime('%Y'), 
                                                          date.strftime('%m')))   
            
            print('Processing {}'.format(date.strftime('%Y-%m')))
            mth = ds.sel(time='{}-{}'.format(date.strftime('%Y'), date.strftime('%m')))
            
            if not os.path.exists(tfile):
                exportTides(mth, tfile)
            else:
                continue
        
        elif scenario == 'era5':
            exp = 'reanalysis'
            
        else:
            full_experiment = raw_data[scenario]['fname'].split('_')[1]
            
            if 'hist' in full_experiment:
                exp = 'historical'
            else:
                exp = 'future'
            
            if year < 2015:
                texp = 'historical'
            else:
                texp = 'future'
                
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
        

def attrib(dataset):
    # add global attribute metadata
    if 'tide' in list(dataset.keys()):
        gvm = 'tide'
    elif 'waterlevel' in list(dataset.keys()):
        gvm = 'waterlevel'
    elif 'surge' in list(dataset.keys()):
        gvm = 'surge'

    if scenario == 'era5':
        experiment = 'reanalysis'
    elif scenario == 'tides':
        experiment = ''
    else:
        experiment = raw_data[scenario]['fname'].split('_')[1]
        
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
                  'institution': 'Deltares', #include partners later on.,
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

raw2nc(year,scenario)

