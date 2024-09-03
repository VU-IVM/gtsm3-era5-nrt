# ---
# Author: N. Aleksandrova
# Contact: natalia.aleksandrova@deltares.nl
# Date created: May 2024
# Remarks:

from datetime import datetime 
import pandas as pd
from global_map import global_map
import warnings
warnings.filterwarnings('ignore')
import xarray as xr
import numpy as np
import sys
import os
import glob
import warnings
import matplotlib as mpl
import matplotlib.ticker as tkr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import cartopy as crt
sys.path.append("..")
from path_dict import path_dict

# PLOTTING inputs
cmap = mpl.colormaps['viridis'].resampled(20)
cmap2 = mpl.colormaps['hot_r']#.resampled(20)
cmap3 = mpl.colormaps['magma_r'].resampled(20)
cmap4 = mpl.colormaps['seismic']#.resampled(20)

# locate timeseries 
dir_data = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/'
dir_wl_cds_raw = os.path.join('/gpfs/work1/0/einf3499/06_model_runs/99_slr_tide_surge_runs_oldproject/')
dir_tide = os.path.join('/gpfs/work1/0/einf3499/03_tides_CDS')
dir_msl_cds = os.path.join(dir_data,'yearly_MSL_from_CDS')

# station ids
stations_gtsm = [18227] # Nigeria
#stations_gtsm = [24113] # Norway
#stations_gtsm = [15349] # Houston

def read_his_file(path):
    ds = xr.open_dataset(path)
    ds = ds.assign_coords({'stations': ds.stations})
    keys = list(ds.keys()) # remove all variables except water level
    keys.remove('waterlevel')
    if 'patm' in keys:
        keys.remove('patm')
    ds = ds.drop(keys)
    return ds

# Open model data
dir_modeloutput = path_dict['modelruns']
path_modeloutput = os.path.join(dir_modeloutput, f'model_input_ERA5_1978','output','gtsm_fine_0000_his.nc')
ds_modelrun = read_his_file(path_modeloutput)
ds_modelrun_sel = ds_modelrun.sel(stations=stations_gtsm,drop=True).squeeze()

## Checking the results for loc 18227


##


### Load raw data from previous (CDS) project
path_modeloutput = os.path.join(dir_wl_cds_raw, f'model_input_ERA5_1979','gtsm_fine_0000_his.nc')
ds_modelrun_cds = read_his_file(path_modeloutput)
ds_modelrun_cds_sel = ds_modelrun_cds.sel(stations=stations_gtsm,drop=True).squeeze()

### Make comparison of single station timeseries
fig, ax = plt.subplots(figsize=(12,7))
ts = ax.plot(ds_modelrun_sel['time'].values,ds_modelrun_sel['waterlevel'].values,'b-',alpha=0.5,label='New model run');
ts3 = ax.plot(ds_modelrun_cds_sel['time'].values,ds_modelrun_cds_sel['waterlevel'].values,'r-',alpha=0.5,label='Old run');
#ts3 = ax.plot(ds_cds_sel['time'].values,ds_cds_sel['waterlevel'].values,'r-',alpha=0.5,label='Data from CDS');
ax.set_ylabel('Water level [mMSL]')
ax.grid(); ax.title.set_text('Timeseries comparison for different runs')

### check difference
print(f'Mean WL 1979/05 in new model runs: {np.nanmean(ds_modelrun["waterlevel"].sel(time=slice("1979-05-01","1979-05-43")).values):.02} m')
print(f'Mean WL 1979/05 in old model runs: {np.nanmean(ds_modelrun_cds["waterlevel"].sel(time=slice("1979-05-01","1979-02-28")).values):.02} m')

### Compare maps of water levels between old and new runs
lims = [-0.5,0.5]

# map of WLs in the new and old runs
fig = plt.figure(figsize=(15,15))
axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 1), share_all=True, axes_pad=1,cbar_location='right',cbar_mode='single',cbar_size='3%',cbar_pad=0.3, label_mode='keep')
ax = global_map(axs[0])
ax = global_map(ax)
bs = ax.scatter(x=ds_modelrun.station_x_coordinate.values,y=ds_modelrun.station_y_coordinate.values,c=ds_modelrun['waterlevel'].sel(time='1979-02-01 00:00:00').values,transform=crt.crs.PlateCarree(),cmap=cmap4,vmin=lims[0],vmax=lims[1]); 
ax.title.set_text(f"Water level - new run - 1 Feb 1979")
ax = global_map(axs[1])
ax = global_map(ax)
bs = ax.scatter(x=ds_modelrun_cds.station_x_coordinate.values,y=ds_modelrun_cds.station_y_coordinate.values,c=ds_modelrun_cds['waterlevel'].sel(time='1979-02-01 00:00:00').values,transform=crt.crs.PlateCarree(),cmap=cmap4,vmin=lims[0],vmax=lims[1]);
ax.title.set_text(f"Water level - old run - 1 Feb 1979")
cbar = ax.cax.colorbar(bs); cbar.set_label('WL [m]',fontsize=14)

# Map of WL difference between new and old runs
fig = plt.figure(figsize=(15,8))
axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(1, 1), share_all=True, axes_pad=1,cbar_location='right',cbar_mode='single',cbar_size='3%',cbar_pad=0.3, label_mode='keep')
ax = global_map(axs[0])
ax = global_map(ax)
bs = ax.scatter(x=ds_modelrun.station_x_coordinate.values,y=ds_modelrun.station_y_coordinate.values,c=ds_modelrun['waterlevel'].sel(time='1979-02-01 00:00:00').values-ds_modelrun_cds['waterlevel'].sel(time='1979-02-01 00:00:00').values,transform=crt.crs.PlateCarree(),cmap=cmap4,vmin=-0.2,vmax=0.2); 
ax.title.set_text(f"Water level difference - new run vs. old run - 1 Feb 1979")
cbar = ax.cax.colorbar(bs); cbar.set_label('WL difference [m]',fontsize=14)


'''
#find station coordinates
lat = 29.3
lon = -91.83

# First, find the index of the grid point nearest a specific lat/lon.   
abslat = np.abs(ds_surge_cds.station_y_coordinate.values-lat)
abslon = np.abs(ds_surge_cds.station_x_coordinate.values-lon)
c = np.maximum(abslon, abslat)

([iloc]) = np.where(c == np.min(c))

# Now I can use that index location to get the values at the x/y diminsion
ds_surge_cds.isel(stations=iloc).load()
'''

### Check tide
path_tide = os.path.join(dir_modeloutput, f'model_input_TIDE_1980_oldSLR','output','gtsm_fine_0000_his.nc')
ds_tide = read_his_file(path_tide)

path_tide_noslr = os.path.join(dir_modeloutput, f'model_input_TIDE_1980_noSLR','output','gtsm_fine_0000_his.nc')
ds_tide_noslr = read_his_file(path_tide)

path_tide = os.path.join('/projects/0/einf3499/06_model_runs/99_slr_tide_surge_runs_oldproject/model_input_tide_1980/','output','gtsm_fine_0000_his.nc')
ds_tide_cds = read_his_file(path_tide)

timestamp = '1978-12-17 00:00:00'

# map of WLs in the new and old runs
fig = plt.figure(figsize=(15,15))
axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 1), share_all=True, axes_pad=1,cbar_location='right',cbar_mode='single',cbar_size='3%',cbar_pad=0.3, label_mode='keep')
ax = global_map(axs[0])
ax = global_map(ax)
bs = ax.scatter(x=ds_tide.station_x_coordinate.values,y=ds_tide.station_y_coordinate.values,c=ds_tide['waterlevel'].sel(time=timestamp).values,transform=crt.crs.PlateCarree(),cmap=cmap4,vmin=lims[0],vmax=lims[1]); 
ax.title.set_text(f"Water level - new run - {timestamp}")
ax = global_map(axs[1])
ax = global_map(ax)
bs = ax.scatter(x=ds_tide_cds.station_x_coordinate.values,y=ds_tide_cds.station_y_coordinate.values,c=ds_tide_cds['waterlevel'].sel(time=timestamp).values,transform=crt.crs.PlateCarree(),cmap=cmap4,vmin=lims[0],vmax=lims[1]);
ax.title.set_text(f"Water level - old run - {timestamp}")
cbar = ax.cax.colorbar(bs); cbar.set_label('WL [m]',fontsize=14)

# Map of tide difference between new and old runs
fig = plt.figure(figsize=(15,8))
axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(1, 1), share_all=True, axes_pad=1,cbar_location='right',cbar_mode='single',cbar_size='3%',cbar_pad=0.3, label_mode='keep')
ax = global_map(axs[0])
ax = global_map(ax)
bs = ax.scatter(x=ds_tide.station_x_coordinate.values,y=ds_tide.station_y_coordinate.values,c=ds_tide['waterlevel'].sel(time=timestamp).values-ds_tide_cds['waterlevel'].sel(time=timestamp).values,transform=crt.crs.PlateCarree(),cmap=cmap4,vmin=-11,vmax=-5); 
ax.title.set_text(f"Tide level difference - new run vs. old run - 1 Feb 1980")
cbar = ax.cax.colorbar(bs); cbar.set_label('Tidel level difference [m]',fontsize=14)

### Make comparison of single station timeseries for tide
fig, ax = plt.subplots(figsize=(12,7))
ts = ax.plot(ds_tide['time'].values,ds_tide['waterlevel'].sel(stations=stations_gtsm).values,'b-',alpha=0.5,label='New model run');
ts3 = ax.plot(ds_tide_cds['time'].values,ds_tide_cds['waterlevel'].sel(stations=stations_gtsm).values,'r-',alpha=0.5,label='Old run');
#ts3 = ax.plot(ds_cds_sel['time'].values,ds_cds_sel['waterlevel'].values,'r-',alpha=0.5,label='Data from CDS');
ax.set_ylabel('Tide level [mMSL]')
ax.grid(); ax.title.set_text('Timeseries comparison for different runs')

### Check another year
path = os.path.join(dir_modeloutput, f'model_input_ERA5_1955','output','gtsm_fine_0000_his.nc')
ds_year = read_his_file(path)


### Open SLR file
ds_slr = xr.open_dataset(os.path.join('/projects/0/einf3499/05_meteo_SLR/','TotalSeaLevel_MapsSROCC_rcp85_Perc50_zero1986to2005_dflow_extrap.nc'))
ds_slr['TotalRelativeSeaLevel'].sel(time='1980').plot()


### Compare post-processed data

""""
### Load CDS results
file_nc = os.path.join(dir_wlts_CDS,f'reanalysis_waterlevel_10min_{year}_{mnth:02}_v1.nc')
ds_cds = xr.open_dataset(file_nc,chunks={'stations': 1000}); 
ds_cds_sel = ds_cds.sel(stations = stations_gtsm,drop=True)

file_nc2 = os.path.join(dir_wlts_CDS,f'reanalysis_waterlevel_10min_2014_02_v1.nc')
ds_cds2 = xr.open_dataset(file_nc2,chunks={'stations': 1000}); 
ds_cds2 = ds_cds2.sel(stations = stations_gtsm,drop=True)
"""" 