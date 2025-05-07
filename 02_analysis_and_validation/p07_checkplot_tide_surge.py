# ---
# Author: N. Aleksandrova
# Contact: natalia.aleksandrova@deltares.nl
# Date created: Oct 2024
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
dir_wl = os.path.join(dir_data,'02_slr_tide_surge_runs')
dir_tide = os.path.join('/gpfs/work1/0/einf3499/03_tides_CDS')

stations_gtsm = 21152 # Vigo, Spain

def read_his_file(path):
    ds = xr.open_dataset(path)
    ds = ds.assign_coords({'stations': ds.stations})
    keys = list(ds.keys()) # remove all variables except water level
    if 'waterlevel' in keys:
        keys.remove('waterlevel')
    if 'tide' in keys:
        keys.remove('tide')        
    if 'patm' in keys:
        keys.remove('patm')
    ds = ds.drop(keys)
    return ds

# Open model data
dir_modeloutput = path_dict['modelruns']
path_modeloutput = os.path.join(dir_modeloutput, f'model_input_ERA5_2023','output','gtsm_fine_0000_his.nc')
ds_era5 = read_his_file(path_modeloutput)
ds_era5_sel = ds_era5.sel(stations=stations_gtsm,time=slice('2023-01-01','2023-01-31'),drop=True).squeeze()

path_tide = os.path.join(dir_modeloutput, f'model_input_TIDE_2023_noSLR','output','gtsm_fine_0000_his.nc')
ds_tide_new = read_his_file(path_tide)
ds_tide_new_sel = ds_tide_new.sel(stations=stations_gtsm,time=slice('2023-01-01','2023-01-31'),drop=True).squeeze()

path_tide = os.path.join(dir_tide, 'future_tide_2023_01_v1.nc')
ds_tide_old = read_his_file(path_tide)
ds_tide_old_sel = ds_tide_old.sel(stations=stations_gtsm,time=slice('2023-01-01','2023-01-31'),drop=True).squeeze()

ds_surge_old = ds_era5_sel['waterlevel'] - ds_tide_old_sel['tide']
ds_surge_new = ds_era5_sel['waterlevel'] - ds_tide_new_sel['waterlevel']

### Make comparison of single station timeseries
fig, ax = plt.subplots(figsize=(12,7))
ts = ax.plot(ds_modelrun_sel['time'].values,ds_modelrun_sel['waterlevel'].values,'b-',alpha=0.5,label='New model run');
ts3 = ax.plot(ds_modelrun_cds_sel['time'].values,ds_modelrun_cds_sel['waterlevel'].values,'r-',alpha=0.5,label='Old run');
#ts3 = ax.plot(ds_cds_sel['time'].values,ds_cds_sel['waterlevel'].values,'r-',alpha=0.5,label='Data from CDS');
ax.set_ylabel('Water level [mMSL]')
ax.grid(); ax.title.set_text('Timeseries comparison for different runs')



