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
    settings = {'yearmin1': 1951,'yearmax1': 2022,
                'yearmin2': 1979,'yearmax2': 2022,
                'yearmin3': 1985,'yearmax3': 2014}
    rps =[2,10,50,100]    
    
    # location of EVA data
    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_eva = os.path.join(dir_postproc,'EVA-GTSM-ERA5')
    dir_eva1 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin1"]}_{settings["yearmax1"]}') #output dir
    dir_eva2 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin2"]}_{settings["yearmax2"]}')
    dir_eva3 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin3"]}_{settings["yearmax3"]}')
    
    #locate .csv files --> TODO: make a function instead? to save repetition
    filenames1 = 'ds_GTSM-ERA5_%s_%s_stations*eva.csv' % (str(settings['yearmin1']),str(settings['yearmax1']))
    dir_data1 = os.path.join(dir_eva1,filenames1)
    file_list1 = glob.glob(dir_data1)
    file_list1.sort()

    filenames2 = 'ds_GTSM-ERA5_%s_%s_stations*eva.csv' % (str(settings['yearmin2']),str(settings['yearmax2']))
    dir_data2 = os.path.join(dir_eva2,filenames2)
    file_list2 = glob.glob(dir_data2)
    file_list2.sort() 

    filenames3 = 'ds_GTSM-ERA5_%s_%s_stations*eva.csv' % (str(settings['yearmin3']),str(settings['yearmax3']))
    dir_data3 = os.path.join(dir_eva3,filenames3)
    file_list3 = glob.glob(dir_data3)
    file_list3.sort()   

    #read and merge eva data in one data frame
    ds_gtsm_eva1 = pd.read_csv(file_list1[0])
    ds_gtsm_eva2 = pd.read_csv(file_list2[0])
    ds_gtsm_eva3 = pd.read_csv(file_list3[0]) #--> for comparison with CDS dataset of RVs
    for ii in range(1,min(len(file_list1),len(file_list2))):
        tmp1 = pd.read_csv(file_list1[ii])
        tmp2 = pd.read_csv(file_list2[ii])
        ds_gtsm_eva1 = pd.concat([ds_gtsm_eva1,tmp1])
        ds_gtsm_eva2 = pd.concat([ds_gtsm_eva2,tmp2])
    del tmp1, tmp2

    for ii in range(1,len(file_list3)):
        tmp = pd.read_csv(file_list3[ii])
        ds_gtsm_eva3 = pd.concat([ds_gtsm_eva3,tmp])
    
    #locate and read files that contain coordinates - TODO: ref to stats files instead of TS files - much smaller files
    file_list_nc = [str(file_list1[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list1))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    coord_x = ds['station_x_coordinate'].values
    coord_y = ds['station_y_coordinate'].values
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        coord_x = np.append(coord_x, ds['station_x_coordinate'].values)
        coord_y = np.append(coord_y, ds['station_y_coordinate'].values)
    
    # Plotting of difference between the two periods
    cmap = mpl.colormaps['viridis'].resampled(20)
    vrange=[0,0.03]    
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
                        s=50,c=ds_gtsm_eva1[str(rps[ii])]-ds_gtsm_eva2[str(rps[ii])],transform=crt.crs.PlateCarree(),
                        cmap=cmap, vmin=vrange[0], vmax=vrange[1])
        cbar = ax.cax.colorbar(bs)
    figname = 'EVA_map_compare_datasets.png' 
    fig.savefig(f'{dir_eva}/{figname}')

    # loading CDS return values
    dir_rv = os.path.join(dir_eva,'EVA_RV_from_CDS')
    ds_rv100 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp100_best-fit_v1.nc'); ds_rv100.close()
    ds_rv50 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp50_best-fit_v1.nc'); ds_rv50.close()
    ds_rv10 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp10_best-fit_v1.nc'); ds_rv10.close()
    ds_rv2 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp2_best-fit_v1.nc'); ds_rv2.close()

    ds_rv = np.array([ds_rv2.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values,
                    ds_rv10.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values,
                     ds_rv50.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values,
                    ds_rv100.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values])
    
    # Plotting difference of present calculation and dataset on CDS
    cmap = mpl.colormaps['magma_r'].resampled(20)
    vrange=[0,0.01]    
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
        bs = ax.scatter(x=coord_x[0:1000],y=coord_y[0:1000],
                        s=50,c=ds_gtsm_eva3[str(rps[ii])]-ds_rv[ii,:],transform=crt.crs.PlateCarree(),
                        cmap=cmap, vmin=vrange[0], vmax=vrange[1])
        cbar = ax.cax.colorbar(bs)
    figname = 'EVA_map_compare_dataset_to_CDS.png' 
    fig.savefig(f'{dir_eva}/{figname}')

