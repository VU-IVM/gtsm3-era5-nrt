# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# Author: Sanne Muis, N. Aleksandrova
# Contact: sanne.muis@deltares.nl
# Date created: July 2023
# Remarks: gtsm_eva

#import modin.pandas as pd
import pandas as pd
from global_map import global_map
import warnings
warnings.filterwarnings('ignore')

#from pathlib import Path
import xarray as xr
import numpy as np
import sys
import os
import glob
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import cartopy as crt

if __name__ == "__main__":   
    settings = {'yearmin': 1951,
                'yearmax': 1978}
    rps =[2,10,50,100]    
    #chunks = 1000 #eva data is saved in chuncks of 1000 stations
    
    # location of EVA data
    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_eva = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin"]}_{settings["yearmax"]}') #output dir
    
    #locate .csv files
    filenames = 'ds_GTSM-ERA5_%s_%s_stations*eva.csv' % (str(settings['yearmin']),str(settings['yearmax']))
    dir_data = os.path.join(dir_eva,filenames)
    file_list = glob.glob(dir_data)
    file_list.sort()

    #read and merge eva data in one data frame
    ds_gtsm_eva = pd.read_csv(file_list[0])
    for ii in range(1,len(file_list)):
        tmp = pd.read_csv(file_list[ii])
        ds_gtsm_eva = pd.concat([ds_gtsm_eva,tmp])
    #eva_data

    #locate and read files that contain coordinates - TODO: ref to stats files instead of TS files - much smaller files
    file_list_nc = [str(file_list[ii]).replace('_eva.csv','.nc') for ii in range(0,len(file_list))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    coord_x = ds['station_x_coordinate'].values
    coord_y = ds['station_y_coordinate'].values
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        coord_x = np.append(coord_x, ds['station_x_coordinate'].values)
        coord_y = np.append(coord_y, ds['station_y_coordinate'].values)
    
    # Plotting
    cmap = mpl.colormaps['viridis'].resampled(20)
    vrange=[0,6]    
    print('Plotting... ')
    fig = plt.figure(figsize=(26,20))
    axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    axs = AxesGrid(fig, 111, axes_class=axes_class,
            nrows_ncols=(2, 2),
            share_all=True,
            axes_pad=1.7,
            cbar_location='right',
            cbar_mode='each',
            cbar_size='3%',
            cbar_pad=1.0,
            label_mode='keep')
    for ii in range(0,4):
        ax = global_map(axs[ii])
        bs = ax.scatter(x=coord_x,y=coord_y,
                        s=50,c=ds_gtsm_eva[str(rps[ii])],transform=crt.crs.PlateCarree(),
                        cmap=cmap, vmin=vrange[0], vmax=vrange[1])
        cbar = ax.cax.colorbar(bs)
    figname = 'EVA_map_%s_%s_stations_ALL.png' % (str(settings['yearmin']),str(settings['yearmax']))
    fig.savefig(f'{dir_eva}/{figname}')
    print('Done!')


