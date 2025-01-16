# ---
# Author: N. Aleksandrova
# Contact: natalia.aleksandrova@deltares.nl
# Date created: May 2024
# Remarks:

import pandas as pd
from global_map import global_map
import warnings
warnings.filterwarnings('ignore')
import xarray as xr
import numpy as np
import sys
import os
from math import sqrt, cos, radians
import glob
import warnings
import matplotlib as mpl
import matplotlib.ticker as tkr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import cartopy as crt
#from pyextremes import EVA
#from pyextremes.plotting import plot_return_values
sys.path.append("..")
from path_dict import path_dict

# function for detrending of timeseries
def detrend(ds: xr.DataArray, plot = False):
  ''' remove annual means and overall mean '''
  ds = ds.assign_coords(year=ds.time.dt.strftime("%Y"))
  ds_new = (ds.groupby("year") - ds.groupby("year").mean("time"))
  ds['sea_level_detrended'] = ds_new['sea_level'] - ds_new['sea_level'].mean()
  if plot == True:
      fig, axs = plt.subplots(nrows=2)
      ds.sea_level.plot.line(x='time',ax=axs[0], add_legend=False)   
      ds.sea_level_detrended.plot.line(x='time',ax=axs[1],add_legend=False)  
  return ds


# PLOTTING inputs
cmap = mpl.colormaps['viridis'].resampled(20)
cmap2 = mpl.colormaps['hot_r']#.resampled(20)
cmap3 = mpl.colormaps['magma_r'].resampled(20)
cmap4 = mpl.colormaps['seismic']#.resampled(20)

# locate timeseries 
dir_data = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/'
dir_wlts = os.path.join(dir_data,'timeseries-GTSM-ERA5-hourly','waterlevel')
dir_wlts_10min = os.path.join(dir_data,'timeseries-GTSM-ERA5-10min','waterlevel')
dir_wlts_old = os.path.join(dir_data,'timeseries-GTSM-ERA5-hourly','waterlevel_old')
dir_wlts2 = os.path.join(dir_data,'timeseries-GTSM-ERA5-hourly-1979-2018','waterlevel')
dir_surgets = os.path.join(dir_data,'timeseries-GTSM-ERA5-hourly','surge')
dir_surgets_old = os.path.join(dir_data,'timeseries-GTSM-ERA5-hourly','surge_old')
dir_surgets2 = os.path.join(dir_data,'timeseries-GTSM-ERA5-hourly-1979-2018','surge')


dir_postproc = path_dict['postproc']
dir_ts = os.path.join(dir_postproc,'timeseries_plots')

''''
dir_wlts = os.path.join(dir_data,'timeseries-GTSM-ERA5-10min','waterlevel')
dir_wlts_old = os.path.join(dir_data,'timeseries-GTSM-ERA5-10min','waterlevel_old')
dir_wlts2 = os.path.join(dir_data,'timeseries-GTSM-ERA5-10min-1979-2018','waterlevel')
dir_surgets = os.path.join(dir_data,'timeseries-GTSM-ERA5-10min','surge')
dir_surgets_old = os.path.join(dir_data,'timeseries-GTSM-ERA5-10min','surge_old')
dir_surgets2 = os.path.join(dir_data,'timeseries-GTSM-ERA5-10min-1979-2018','surge')
''''

# station ids
#stations_gtsm = [18227] # Nigeria
#stations_gtsm = [24113] # Norway
#stations_gtsm = [15349] # Houston

stations_gtsm = [15349,24113,18227,17584,38671,41030,38637,40564,40353,39598,17432,15160,15759]
'''
file_nc = os.path.join(dir_surgets,f'reanalysis_surge_hourly_1978_01_v1.nc')        
ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); 

#find station coordinates
lat = 4.52
lon = 8.4

# First, find the index of the grid point nearest a specific lat/lon.   
abslat = np.abs(ds.station_y_coordinate.values-lat)
abslon = np.abs(ds.station_x_coordinate.values-lon)
c = np.maximum(abslon, abslat)

([iloc]) = np.where(c == np.min(c))

station = ds.stations.values[iloc[0]]
'''

year_start = 2023
year_end = 2025

#load GTSM timeseries for selected stations
for year in range(year_start,year_end):
    print(f'loading {year}')
    for mnth in range(1,13):
        print(f'loading {year}, month {mnth}')
        if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
            file_nc = os.path.join(dir_surgets,f'reanalysis_surge_hourly_{year}_{mnth:02d}_v1.nc')
            #file_nc = os.path.join(dir_surgets,f'reanalysis_surge_10min_{year}_{mnth:02d}_v1.nc')
            ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
            #ds = ds.drop('station_name')
            #file_nc_old = os.path.join(dir_surgets_old,f'reanalysis_surge_hourly_{year}_{mnth:02d}_v1.nc')
            #file_nc_old = os.path.join(dir_surgets_old,f'reanalysis_surge_10min_{year}_{mnth:02d}_v1.nc')
            #ds_old = xr.open_dataset(file_nc_old,chunks={'stations': 1000}); 
            #ds_old=ds_old.sel(stations=stations_gtsm,drop=True)
            #ds_old.load()
            
            if "station_x_coordinate" in list(ds.data_vars):
                ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))
            #if "station_x_coordinate" in list(ds_old.data_vars):
                #ds_old = ds_old.set_coords(("station_x_coordinate", "station_y_coordinate"))
        else:
            #file_nc = os.path.join(dir_surgets2,f'reanalysis_surge_10min_{year}_{mnth:02d}_v1.nc')
            file_nc = os.path.join(dir_surgets2,f'reanalysis_surge_hourly_{year}_{mnth:02d}_v1.nc')
            ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
        
        ds=ds.sel(stations=stations_gtsm,drop=True)
        ds.load()
        
        if ((year == year_start) & (mnth == 1)):
            ds_gtsm = ds
            #ds_gtsm_old = ds_old
        else:
            ds_gtsm = xr.concat([ds_gtsm,ds],dim="time")
            
            #if ((year < 1980) | (year > 2018)):  
            #    ds_gtsm_old = xr.concat([ds_gtsm_old,ds_old],dim="time")
            #    del ds_old
        del ds

# open timeseries selection and detrend
#ds_gtsm['sea_level'] = ds_gtsm['waterlevel']
#ds_gtsm = ds_gtsm.drop(['waterlevel'])
#ds_gtsm = detrend(ds_gtsm)
#ds_gtsm = ds_gtsm.drop(['sea_level'])
#ds_gtsm = ds_gtsm.chunk({"time": -1, "stations": "auto"})
#ds_gtsm.load()
#ds_gtsm=ds_gtsm.set_coords(("station_x_coordinate", "station_y_coordinate"))

# make overview plots per location
for ss in range(5,len(stations_gtsm)):

    print('processing station ',ss,' out of ', len(stations_gtsm))

    ts_gtsm = ds_gtsm.surge.sel(stations=stations_gtsm[ss])#.sel(time=slice('01-01-1950','31-12-2022'))
    #ts_gtsm_old = ds_gtsm_old.surge.sel(stations=stations_gtsm[ss]).sel(time=slice('01-01-1950','31-12-2022'))
    

    # plot
    fig = plt.figure(figsize=(20,20))
    ax0 = plt.subplot2grid((2, 2), (0, 0), colspan=1, projection=crt.crs.Robinson())
    ax1 = plt.subplot2grid((2, 2), (1, 0), colspan=2, rowspan=1)

    #    plot location on the map
    ax0 = global_map(ax0)
    bs = ax0.scatter(x=ts_gtsm.station_x_coordinate.values,y=ts_gtsm.station_y_coordinate.values,marker ='.',s=200,transform=crt.crs.PlateCarree(),facecolors='none',edgecolors='red',linewidth=3); 
    ax0.title.set_text(f"Location {int(stations_gtsm[ss])}")

    # plot timeseries
    ts = ax1.plot(ts_gtsm['time'].values,ts_gtsm.values,'b-',alpha=0.5,label='GTSM-ERA5');
    #ts = ax1.plot(ts_gtsm['time'].sel(time=slice('1978-01-01','2000-01-01')).values,ts_gtsm.sel(time=slice('1978-01-01','2000-01-01')).values,'b-',alpha=0.5,label='GTSM-ERA5');
    #ts2 = ax1.plot(ts_gtsm_old['time'].values,ts_gtsm_old.values,'r-',alpha=0.5,label='GTSM-ERA5');
    ax1.set_ylabel('Surge level [m]')
    ax1.grid(); ax1.title.set_text('Full timeseries')


    figname = f'TS_station_{str(stations_gtsm[ss]).zfill(5)}_surge_{year_start}_{year_end}.png' 
    fig.savefig(f'{dir_ts}/{figname}')



### test plot
st_id = 18227
file_nc = os.path.join(dir_surgets,f'reanalysis_surge_hourly_1978_01_v1.nc')
ds_surge_1978 = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
ds_surge_1978 = ds_surge_1978.sel(stations=st_id)

file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_1978_01_v1.nc')
ds_wl_1978 = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
ds_wl_1978 = ds_wl_1978.sel(stations=st_id)

file_nc = os.path.join('/gpfs/work1/0/einf3499/03_tides_CDS',f'historical_tide_1978_01_v1.nc')
ds_tide_1978 = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
ds_tide_1978 = ds_tide_1978.sel(stations=st_id)


file_nc = os.path.join(dir_surgets,f'reanalysis_surge_hourly_1979_01_v1.nc')
ds_surge_1979 = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
ds_surge_1979 = ds_surge_1979.sel(stations=st_id)

file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_1979_01_v1.nc')
ds_wl_1979 = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
ds_wl_1979 = ds_wl_1979.sel(stations=st_id)

file_nc = os.path.join('/gpfs/work1/0/einf3499/03_tides_CDS',f'historical_tide_1979_01_v1.nc')
ds_tide_1979 = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
ds_tide_1979 = ds_tide_1979.sel(stations=st_id)

fig, ax = plt.subplots(figsize=(15,8))
ax.plot(ds_surge_1979.time,ds_surge_1979.surge)
ax.plot(ds_wl_1979.time,ds_wl_1979.waterlevel)
ax.plot(ds_tide_1979.time,ds_tide_1979.tide)
ax.title.set_text(f'station {st_id} - timeseries of waterlevel, surge and tide')
figname = f'TS_station_{str(st_id).zfill(5)}_1979.png' 
fig.savefig(f'{dir_ts}/{figname}')

fig, ax = plt.subplots(figsize=(15,8))
ax.plot(ds_surge_1978.time,ds_surge_1978.surge)
ax.plot(ds_wl_1978.time,ds_wl_1978.waterlevel)
ax.plot(ds_tide_1978.time,ds_tide_1978.tide)
ax.title.set_text(f'station {st_id} - timeseries of waterlevel, surge and tide')
figname = f'TS_station_{str(st_id).zfill(5)}_1978.png' 
fig.savefig(f'{dir_ts}/{figname}')

