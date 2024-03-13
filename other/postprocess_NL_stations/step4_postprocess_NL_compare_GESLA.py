# -*- coding: utf-8 -*-
"""
Created on Mon Aug 01 2023

@author: n-aleksandrova
"""

import xarray as xr
import sys
import os
import matplotlib.pyplot as plt
year = int(sys.argv[1])
print(f'Year: {year}')

from gesla import GeslaDataset

sys.path.append("..")
from path_dict import path_dict
dir_postproc = path_dict['postproc']
dir_gesla = path_dict['gesla']
dir_data = os.path.join(dir_postproc,'NL_stations','timeseries') 
dir_output = os.path.join(dir_postproc,'NL_stations','comparison_GESLA') 

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

station_list = ['NWS_NO_TS_MO_Vlissingen','NWS_NO_TS_MO_HoekVanHolland','NWS_NO_TS_MO_Ijmuiden','NWS_NO_TS_MO_DenHelder','NWS_NO_TS_MO_Harlingen','NWS_NO_TS_MO_Delfzijl']
gesla_list = ['vlissingen-vlissgn-nld-rws','hoek_van_holland-hoekvhld-nld-rws','ijmuiden_noordersluis-ijmdndss-nld-rws','den_helder-denhdr-nld-rws','harlingen-harlgn-nld-rws','delfzijl-delfzl-nld-rws']

# get timeseries for water levels
file_nc = os.path.join(f'{dir_data}/era5_reanalysis_waterlevel_{year}_NLstations.nc')
ds_wl = xr.open_dataset(file_nc); 
    
# get timeseries of surge levels
#file_nc = os.path.join(f'{dir_data}/era5_reanalysis_surge_{year}_NLstations.nc')
#ds_surge = xr.open_dataset(file_nc); 

# get GESLA data
meta_file = os.path.join(dir_gesla,"GESLA3_ALL.csv")
data_path = os.path.join(dir_gesla,"GESLA3.0_ALL//")
filenames = os.listdir(data_path)
obj_gesla = GeslaDataset(meta_file=meta_file, data_path=data_path)
gesla_data = obj_gesla.files_to_xarray(gesla_list)
print(gesla_data)

# plot total water level GTSM and observations
for ii in range(6):
    fig = plt.gcf()
    fig.set_size_inches(15, 7)

    ds_wl_sel = ds_wl.sel(time=slice(f'{year}-01-10',f'{year}-01-25'),stations=station_list[ii])
    plt.plot(ds_wl_sel.time,ds_wl_sel.waterlevel,label='GTSM-ERA5-E')

    gesla_sel = gesla_data.sel(date_time=slice(f'{year}-01-10',f'{year}-01-25'),station=ii)
    plt.plot(gesla_sel.date_time,gesla_sel.sea_level,'ro--',label='Observations')

    plt.xlim(ds_wl_sel.time[0],ds_wl_sel.time[-1])
    plt.ylabel('Water level [m]')
    plt.xticks(rotation = 45) 
    plt.legend(loc='lower right')
    plt.title(f'Water levels at {station_list[ii]}')
    fig.subplots_adjust(bottom=0.15,top=0.9)
    plt.savefig(f'{dir_output}/era5_reanalysis_waterlevel_{year}_{station_list[ii][13:]}_withObs_timeseries.png',bbox_inches='tight')
    plt.clf()

# compare monthly means
ds_wl_mnth = ds_wl.groupby('time.month').mean()
gesla_year = gesla_data.sel(date_time=slice(f'{year}-01-01',f'{year}-12-31'))
gesla_mnth = gesla_year.groupby('date_time.month').mean()

print(ds_wl_mnth)
print(gesla_mnth)

for ii in range(6):
    fig = plt.gcf()
    fig.set_size_inches(15, 7)

    ds_wl_sel = ds_wl_mnth.sel(stations=station_list[ii])
    plt.plot(ds_wl_sel.month,ds_wl_sel.waterlevel,label='GTSM-ERA5-E')

    gesla_sel = gesla_mnth.sel(station=ii)
    plt.plot(gesla_sel.month,gesla_sel.sea_level,'ro--',label='Observations')

    plt.ylabel('Water level [m]')
    plt.xticks(rotation = 45) 
    plt.legend(loc='lower right')
    plt.title(f'Water level monthly means in {year} at {station_list[ii]}')
    fig.subplots_adjust(bottom=0.15,top=0.9)
    plt.savefig(f'{dir_output}/era5_reanalysis_waterlevel_{year}_{station_list[ii][13:]}_withObs_monthly.png',bbox_inches='tight')
    plt.clf()