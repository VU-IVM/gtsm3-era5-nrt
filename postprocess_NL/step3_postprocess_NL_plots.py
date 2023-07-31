# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 2023

@author: n-aleksandrova,veenstra
"""

import xarray as xr
import sys
import os
import matplotlib.pyplot as plt
year = int(sys.argv[1])
print(year)

sys.path.append("..")
from path_dict import path_dict
dir_postproc = path_dict['postproc']
dir_data = os.path.join(dir_postproc,'NL_stations','timeseries') 

dir_output = os.path.join(dir_postproc,'NL_stations','plots') 

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

station_list = ['NWS_NO_TS_MO_Vlissingen','NWS_NO_TS_MO_HoekVanHolland','NWS_NO_TS_MO_Ijmuiden','NWS_NO_TS_MO_DenHelder','NWS_NO_TS_MO_Harlingen','NWS_NO_TS_MO_Delfzijl']

# get timeseries for water levels
file_nc = os.path.join(f'{dir_data}/era5_reanalysis_waterlevel_{year}_NLstations.nc')
ds_wl = xr.open_dataset(file_nc); 
    
# get timeseries of surge levels
#file_nc = os.path.join(f'{dir_data}/era5_reanalysis_surge_{year}_NLstations.nc')
#ds_surge = xr.open_dataset(file_nc); 

plt.figure(figsize=(20,6))

for ii in range(6):
    ds_wl_sel = ds_wl.sel(time=slice('1951-01-01','1951-01-31'),stations=station_list[ii])
    plt.plot(ds_wl_sel.time,ds_wl_sel.waterlevel)

plt.ylabel('Water level [m]')
plt.xticks(rotation = 45) 

plt.savefig(f'{dir_output}/era5_reanalysis_waterlevel_{year}_TS_plot.png',bbox_inches='tight')
    

