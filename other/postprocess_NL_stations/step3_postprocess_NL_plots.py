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
file_nc = os.path.join(f'{dir_data}/era5_reanalysis_surge_{year}_NLstations.nc')
ds_surge = xr.open_dataset(file_nc); 

# plot total water level and surge per station
for ii in range(6):
    fig = plt.gcf()
    fig.set_size_inches(15, 7)

    ds_wl_sel = ds_wl.sel(time=slice(f'{year}-01-15',f'{year}-01-31'),stations=station_list[ii])
    plt.plot(ds_wl_sel.time,ds_wl_sel.waterlevel,label='Total water level')

    ds_surge_sel = ds_surge.sel(time=slice(f'{year}-01-15',f'{year}-01-31'),stations=station_list[ii])
    plt.plot(ds_surge_sel.time,ds_surge_sel.surge,label='Surge')

    plt.xlim(ds_wl_sel.time[0],ds_wl_sel.time[-1])
    plt.ylabel('Water level [m]')
    plt.xticks(rotation = 45) 
    plt.legend(loc='lower right')
    plt.title(f'Total water level at {station_list[ii]}')
    fig.subplots_adjust(bottom=0.15,top=0.9)
    plt.savefig(f'{dir_output}/era5_reanalysis_levels_{year}_TS_{station_list[ii]}.png',bbox_inches='tight')
    plt.clf()
