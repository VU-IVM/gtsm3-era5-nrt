#!/usr/bin/env python
import xarray as xr
import os
import matplotlib.pyplot as plt
import sys

sys.path.append("..")
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from path_dict import path_dict
dir_output = path_dict['meteo_fm']

for year in range(1951,1952):

    for par in ['msl','u10','v10']:

        file_out = os.path.join(dir_output,f'ERA5_CDS_atm_{par}_{year-1}-12-15_{year+1}-01-01.nc')

        with xr.open_dataset(file_out) as data_xr_check:

            for varkey in data_xr_check.data_vars:

                for tm in [0,1,2,25,1000]:

                    print(f'plotting {varkey}')

                    fig,ax1 = plt.subplots()

                    data_xr_check[varkey].isel(time=tm).plot(ax=ax1)

                    fig.savefig(os.path.join(dir_output,'checkplots',f'Plot_{varkey}_{year}_timestamp_{tm}'))
