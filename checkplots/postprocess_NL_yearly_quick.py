# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 16:34:27 2023

@author: veenstra
"""

import xarray as xr
#import matplotlib.pyplot as plt
import sys
year = int(sys.argv[1])
print(year)

#grep -ir 'vlis' /gpfs/work1/0/einf3499/model_runs_extended/slr_tide_surge_runs/model_input_ERA5_1970/selected_output_new_unique.xyn
station_list = ['NWS_NO_TS_MO_Vlissingen','NWS_NO_TS_MO_HoekVanHolland','NWS_NO_TS_MO_Ijmuiden','NWS_NO_TS_MO_DenHelder','NWS_NO_TS_MO_Harlingen','NWS_NO_TS_MO_Delfzijl']

file_nc = f"/gpfs/work1/0/einf3499/model_runs_extended/slr_tide_surge_runs/model_input_ERA5_{year}/output/gtsm_fine_0000_his.nc"
ds = xr.open_dataset(file_nc,chunks={'stations':1000})

#set station_name variable/coordinate for stations dimensions as its index
ds['station_name'] = ds['station_name'].load().str.decode('utf-8',errors='ignore').str.strip() #.load() is essential to convert not only first letter of string.
ds = ds.set_index({'stations':'station_name'})

#drop duplicate indices (stations/crs/gs), this avoids "InvalidIndexError: Reindexing only valid with uniquely valued Index objects"
duplicated_keepfirst = ds['stations'].to_series().duplicated(keep='first')
if duplicated_keepfirst.sum()>0:
    print(f'dropping {duplicated_keepfirst.sum()} duplicate "station_name" labels to avoid InvalidIndexError')
    ds = ds[{'stations':~duplicated_keepfirst}]

data_xr_sel = ds.sel(time=str(year),stations=station_list)
data_xr_yearmean = data_xr_sel.waterlevel.mean(dim='time')
data_xr_yearmean_pd = data_xr_yearmean.to_pandas().round(4)
data_xr_yearmean_pd.index.name = year
print(data_xr_yearmean_pd)
data_xr_yearmean_pd.to_csv(f'yearmean_{year}_sixstations.csv')
