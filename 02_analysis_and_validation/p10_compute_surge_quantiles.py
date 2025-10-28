# ---
# Author: N. Aleksandrova
# Date created: May 2025
# Remarks: This script is for producing various global and location-specific plots of water levels based on the GTSM-ERA5 model runs. 

import warnings
warnings.filterwarnings('ignore')
import xarray as xr
import numpy as np
import sys
import os
import glob
import warnings

# choose which plots to make
make_global_plots = 1
make_local_plots = 1

# location of data
sys.path.append("..")
from path_dict import path_dict
dir_postproc = path_dict['postproc']

dir_surge = os.path.join(dir_postproc, 'timeseries-GTSM-ERA5-hourly', 'surge')
dir_surge2 = os.path.join(dir_postproc, 'timeseries-GTSM-ERA5-hourly-1979-2018', 'surge')
file_list_nc = glob.glob(os.path.join(dir_surge,'*.nc'))
file_list_nc = file_list_nc + glob.glob(os.path.join(dir_surge2,'*.nc'))
file_list_nc.sort()

# get station coordinates
ds = xr.open_dataset(file_list_nc[0]); 
coord_x = ds['station_x_coordinate'].values
coord_y = ds['station_y_coordinate'].values
stations = ds['stations'].values
ds.close()

ds = xr.open_mfdataset(file_list_nc, coords='minimal', compat='override', join='exact')
ds_ext = ds.sel(time=slice('01-01-1950','31-12-1979'))
for ii in range(0,44000,1000): 
    print(f"Processing stations {ii} to {ii+999}")
    ds_temp = ds_ext.sel(stations=slice(ii,ii+999), drop=True)
    ds_temp.load()
    ds_q_sel = ds_temp.surge.quantile([0.90,0.95,0.99], dim=('time'), skipna=False)
    del ds_temp
    ds_q_sel.to_netcdf(os.path.join(dir_postproc,'EVA-GTSM-ERA5','surge_stats',f'GTSM_ERA5_quantiles_1950_1979_stations_{ii:05}_{(ii+999):05}.nc'))
    del ds_q_sel

ds_ori = ds.sel(time=slice('01-01-1990','31-12-2019'))
for ii in range(0,44000,1000): 
    print(f"Processing stations {ii} to {ii+999}")
    ds_temp = ds_ori.sel(stations=slice(ii,ii+999), drop=True)
    ds_temp.load()
    ds_q_sel = ds_temp.surge.quantile([0.90,0.95,0.99], dim=('time'), skipna=False)
    del ds_temp
    ds_q_sel.to_netcdf(os.path.join(dir_postproc,'EVA-GTSM-ERA5','surge_stats',f'GTSM_ERA5_quantiles_1990_2019_stations_{ii:05}_{(ii+999):05}.nc'))
    del ds_q_sel