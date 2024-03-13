# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 16:34:27 2023

@author: veenstra,n-aleksandrova
"""

import xarray as xr
#import matplotlib.pyplot as plt
import sys
import os
year = int(sys.argv[1])
print(year)

sys.path.append("..")
from path_dict import path_dict
dir_runs = path_dict['modelruns']
dir_postproc = path_dict['postproc']
dir_output = os.path.join(dir_postproc,'NL_stations','stats') 

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

station_list = ['NWS_NO_TS_MO_Vlissingen','NWS_NO_TS_MO_HoekVanHolland','NWS_NO_TS_MO_Ijmuiden','NWS_NO_TS_MO_DenHelder','NWS_NO_TS_MO_Harlingen','NWS_NO_TS_MO_Delfzijl']

file_nc = os.path.join(dir_runs,f'model_input_ERA5_{year}/output/gtsm_fine_0000_his.nc')
ds = xr.open_dataset(file_nc,chunks={'stations':1000})

#set station_name variable/coordinate for stations dimensions as its index
ds['station_name'] = ds['station_name'].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.
ds = ds.set_index({'stations':'station_name'})

#drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
duplicated_keepfirst = ds['stations'].to_series().duplicated(keep='first')
if duplicated_keepfirst.sum()>0:
    print(f'dropping {duplicated_keepfirst.sum()} duplicate "station_name" labels to avoid InvalidIndexError')
    ds = ds[{'stations':~duplicated_keepfirst}]

# select specific stations
data_xr_sel = ds.sel(time=str(year),stations=station_list)

# create a dictionary of stats and loop over three parameters to write csv files
data = {}
tmp = data_xr_sel.waterlevel.mean(dim='time')
data['mean'] = tmp
tmp = data_xr_sel.waterlevel.min(dim='time')
data['min'] = tmp
tmp = data_xr_sel.waterlevel.max(dim='time')
data['max'] = tmp

for key, value in data.items():
    data_pd = xr.DataArray.to_pandas(value)
    data_pd = data_pd.round(4)
    data_pd.index.name = year
    print(data_pd)
    data_pd.to_csv(f'{dir_output}/Stations_NL_{year}_{key}.csv')