# ---
# Author: N. Aleksandrova
<<<<<<< HEAD
=======
# Contact: natalia.aleksandrova@deltares.nl
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
# Date created: July 2023
# Remarks: This script is for producing various global and location-specific plots of water levels based on the GTSM-ERA5 model runs. 

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
from pyextremes import EVA
from pyextremes.plotting import plot_return_values

# choose which plots to make
make_global_plots = 0
make_local_plots = 1

<<<<<<< HEAD
# function for detrending timeseries
=======
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
def detrend(ds: xr.DataArray, plot = False):
  ''' remove annual means and overall mean '''
  ds = ds.assign_coords(year=ds.time.dt.strftime("%Y"))
  ds_new = (ds.groupby("year") - ds.groupby("year").mean("time"))
  ds['sea_level_detrended'] = ds_new['sea_level'] - ds_new['sea_level'].mean()
  if plot == True:
      fig, axs = plt.subplots(nrows=2)
      ds.sea_level.plot.line(x='time',ax=axs[0], add_legend=False)   
      ds.sea_level_detrended.plot.line(x='time',ax=axs[1],add_legend=False)  
  return ds

if __name__ == "__main__":   
    # EVA was made available for three periods:
<<<<<<< HEAD
    settings = {'yearmin0': 1979,'yearmax0': 2018,'mode0':'1hr',
                'yearmin1': 1950,'yearmax1': 2022,'mode1':'1hr'}
=======
    settings = {'yearmin1': 1985,'yearmax1': 2014,'mode1':'1hr',
                'yearmin2': 1979,'yearmax2': 2018,'mode2':'1hr',
                'yearmin3': 1950,'yearmax3': 2022,'mode3':'1hr'}
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
    rps =[1,10,50,100]    
    
    # location of EVA data
    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_eva = os.path.join(dir_postproc,'EVA-GTSM-ERA5')
<<<<<<< HEAD
    
    #specify directories where timeseries data is stored
    dir_eva0 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin0"]}_{settings["yearmax0"]}_{settings["mode0"]}_v2')
    dir_eva1 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin1"]}_{settings["yearmax1"]}_{settings["mode1"]}_v2')
    
    #locate .csv files
    filenames0 = 'ds_GTSM-ERA5_%s_%s_stations*eva.csv' % (str(settings['yearmin0']),str(settings['yearmax0']))
    dir_data0 = os.path.join(dir_eva0,filenames0)
    file_list0 = glob.glob(dir_data0)
    file_list0.sort() 

    filenames1 = 'ds_GTSM-ERA5_%s_%s_stations*eva.csv' % (str(settings['yearmin1']),str(settings['yearmax1']))
    dir_data1 = os.path.join(dir_eva1,filenames1)
    file_list1 = glob.glob(dir_data1)
    file_list1.sort()   

    #read and merge eva data in one data frame
    ds_gtsm_eva0 = pd.read_csv(file_list0[0])
    ds_gtsm_eva1 = pd.read_csv(file_list1[0]) 
    for ii in range(1,min(len(file_list0),len(file_list1))):
        tmp2 = pd.read_csv(file_list0[ii])
        ds_gtsm_eva0 = pd.concat([ds_gtsm_eva0,tmp2])
        tmp = pd.read_csv(file_list1[ii])
        ds_gtsm_eva1 = pd.concat([ds_gtsm_eva1,tmp])
    del tmp, tmp2
        
    #locate and read files that contain coordinates 
    file_list_nc = [str(file_list1[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list1))]
=======
    dir_eva1 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin1"]}_{settings["yearmax1"]}_{settings["mode1"]}_v2') #output dir
    dir_eva2 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin2"]}_{settings["yearmax2"]}_{settings["mode2"]}_v2')
    dir_eva3 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin3"]}_{settings["yearmax3"]}_{settings["mode3"]}_v2')
    
    #locate .csv files
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
    ds_gtsm_eva3 = pd.read_csv(file_list3[0]) 
    for ii in range(1,min(len(file_list1),len(file_list2))):
        tmp1 = pd.read_csv(file_list1[ii])
        tmp2 = pd.read_csv(file_list2[ii])
        ds_gtsm_eva1 = pd.concat([ds_gtsm_eva1,tmp1])
        ds_gtsm_eva2 = pd.concat([ds_gtsm_eva2,tmp2])
    del tmp1, tmp2

    for ii in range(1,len(file_list3)):
        tmp = pd.read_csv(file_list3[ii])
        ds_gtsm_eva3 = pd.concat([ds_gtsm_eva3,tmp])
    
    #locate and read files that contain coordinates 
    file_list_nc = [str(file_list3[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list3))]
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    coord_x = ds['station_x_coordinate'].values
    coord_y = ds['station_y_coordinate'].values
    stations = ds['stations'].values
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        coord_x = np.append(coord_x, ds['station_x_coordinate'].values)
        coord_y = np.append(coord_y, ds['station_y_coordinate'].values)
        stations = np.append(stations, ds['stations'].values)

    # read quantiles
<<<<<<< HEAD
    file_list_nc = [str(file_list0[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list0))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    quan90_0 = ds['sea_level_detrended'].values[5,:] # 5th percentile is 90% in this dataset, check if correct
    quan95_0 = ds['sea_level_detrended'].values[6,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        quan90_0 = np.append(quan90_0,ds['sea_level_detrended'].values[5,:])
        quan95_0 = np.append(quan95_0,ds['sea_level_detrended'].values[6,:])
    del ds, file_list_nc

=======
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
    file_list_nc = [str(file_list1[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list1))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    quan90_1 = ds['sea_level_detrended'].values[5,:]
    quan95_1 = ds['sea_level_detrended'].values[6,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        quan90_1 = np.append(quan90_1,ds['sea_level_detrended'].values[5,:])
        quan95_1 = np.append(quan95_1,ds['sea_level_detrended'].values[6,:])
<<<<<<< HEAD
=======
    del ds, file_list_nc
    
    file_list_nc = [str(file_list2[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list2))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    quan90_2 = ds['sea_level_detrended'].values[5,:]
    quan95_2 = ds['sea_level_detrended'].values[6,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        quan90_2 = np.append(quan90_2,ds['sea_level_detrended'].values[5,:])
        quan95_2 = np.append(quan95_2,ds['sea_level_detrended'].values[6,:])
    del ds, file_list_nc

    file_list_nc = [str(file_list3[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list3))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    quan90_3 = ds['sea_level_detrended'].values[5,:]
    quan95_3 = ds['sea_level_detrended'].values[6,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        quan90_3 = np.append(quan90_3,ds['sea_level_detrended'].values[5,:])
        quan95_3 = np.append(quan95_3,ds['sea_level_detrended'].values[6,:])
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
    del ds, file_list_nc    

    # PLOTTING
    print('Plotting... ')
    cmap = mpl.colormaps['viridis'].resampled(20)
    cmap2 = mpl.colormaps['hot_r']#.resampled(20)
    cmap3 = mpl.colormaps['magma_r'].resampled(20)
    cmap4 = mpl.colormaps['seismic']#.resampled(20)
    
    if make_global_plots:
        
<<<<<<< HEAD
        # # Plotting percentiles - absolute and difference between periods
        quan95 = np.vstack((quan95_0,quan95_1))
        quan90 = np.vstack((quan90_0,quan90_1))
        vrange1=[0,3]; vrange2=[-0.03,0.03]; 
          
        
=======
        # # Plotting percentiles - periods 1,2,3 and difference
        quan95 = np.vstack((quan95_1,quan95_2,quan95_3))
        quan90 = np.vstack((quan90_1,quan90_2,quan90_3))
        vrange1=[0,3]; vrange2=[-0.03,0.03];   
    # fig = plt.figure(figsize=(26,20))
    # axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    # axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 2), share_all=True, axes_pad=1.7,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=1.0, label_mode='keep')
    # for ii in range(0,3): 
    #     ax = global_map(axs[ii])
    #     bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[ii,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); cbar = ax.cax.colorbar(bs); 
    #     ax.title.set_text(f"95th percentile of total water levels \n Model period: {settings[f'yearmin{ii+1}']}-{settings[f'yearmax{ii+1}']}")
    # ax = global_map(axs[3])
    # bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[2,:]-quan95[1,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); cbar = ax.cax.colorbar(bs); 
    # ax.title.set_text(f"95th percentile - abs. difference, {settings['yearmin3']}-{settings['yearmax3']} vs. {settings['yearmin2']}-{settings['yearmax2']}")
    # figname = 'EVA_map_95percentile_comparison_between_periods_v4.png' 
    # fig.savefig(f'{dir_eva}/{figname}')

        # ---- Smaller plot with just two panes
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
        fig = plt.figure(figsize=(12,8))
        axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
        axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 1), share_all=True, axes_pad=1,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=0.3, label_mode='keep')
        ax = global_map(axs[0])
<<<<<<< HEAD
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[1,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
        cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=14)
        ax.set_title(f"95th percentile of still water levels \n based on {settings[f'yearmin1']}-{settings[f'yearmax1']} dataset",fontsize=16)
        ax = global_map(axs[1])
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[1,:]-quan95[0,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); 
        cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=14)
        ax.set_title(f"Difference in 95th percentile values \n between {settings['yearmin1']}-{settings['yearmax1']} vs. {settings['yearmin0']}-{settings['yearmax0']} dataset",fontsize=16)
=======
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[2,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
        cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=14)
        ax.set_title(f"95th percentile of still water levels \n based on {settings[f'yearmin3']}-{settings[f'yearmax3']} dataset",fontsize=16)
        ax = global_map(axs[1])
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[2,:]-quan95[1,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); 
        cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=14)
        ax.set_title(f"Difference in 95th percentile values \n between {settings['yearmin3']}-{settings['yearmax3']} vs. {settings['yearmin2']}-{settings['yearmax2']} dataset",fontsize=16)
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
        figname = 'EVA_map_95percentile_comparison_between_periods_v5.png'  
        fig.savefig(f'{dir_eva}/{figname}')


<<<<<<< HEAD
        # plots of extreme return periods and difference between periods
        for rp in [1,10,100]:
            rp_all = np.vstack((np.array(ds_gtsm_eva0[f'{str(rp)}_bf']),np.array(ds_gtsm_eva1[f'{str(rp)}_bf'])))
=======
    # #Plotting return values
    # for rp in [1,10,100]:
    #     rp_all = np.vstack((np.array(ds_gtsm_eva1[f'{str(rp)}_bf']),np.array(ds_gtsm_eva2[f'{str(rp)}_bf']),np.array(ds_gtsm_eva3[f'{str(rp)}_bf'])))
    #     vrange1=[0,3]; vrange2=[-0.5,0.5];   
    #     fig = plt.figure(figsize=(26,20))
    #     axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    #     axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 2), share_all=True, axes_pad=1.7,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=1.0, label_mode='keep')

    #     ax = global_map(axs[0])
    #     bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); cbar = ax.cax.colorbar(bs); 
    #     ax.title.set_text(f"Total water level RP{rp}, period {settings[f'yearmin{2}']}-{settings[f'yearmax{2}']}")

    #     ax = global_map(axs[1])
    #     bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); cbar = ax.cax.colorbar(bs); 
    #     ax.title.set_text(f"Total water level RP{rp}, period {settings[f'yearmin{3}']}-{settings[f'yearmax{3}']}")        

    #     ax = global_map(axs[2])
    #     bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:]-rp_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); cbar = ax.cax.colorbar(bs); 
    #     ax.title.set_text(f"Total water level RP{rp}, difference {settings['yearmin3']}-{settings['yearmax3']} vs. {settings['yearmin1']}-{settings['yearmax1']}")

    #     ax = global_map(axs[3])
    #     bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:]-rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); cbar = ax.cax.colorbar(bs); 
    #     ax.title.set_text(f"Total water level RP{rp}, diff. {settings['yearmin3']}-{settings['yearmax3']} vs. {settings['yearmin2']}-{settings['yearmax2']}")        
    #     figname = f'EVA_map_RP{rp}_comparison_between_periods_v4.png' 
    #     fig.savefig(f'{dir_eva}/{figname}')

        # smaller plot with just two panels of extreme return periods
        for rp in [1,10,100]:
            rp_all = np.vstack((np.array(ds_gtsm_eva1[f'{str(rp)}_bf']),np.array(ds_gtsm_eva2[f'{str(rp)}_bf']),np.array(ds_gtsm_eva3[f'{str(rp)}_bf'])))
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            vrange1=[0,5]; vrange2=[-0.5,0.5];  
            fig = plt.figure(figsize=(12,8))
            axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
            axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 1), share_all=True, axes_pad=1.0,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=0.3, label_mode='keep')

            ax = global_map(axs[0])
<<<<<<< HEAD
            bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap3, vmin=vrange1[0], vmax=vrange1[1]); 
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=14)
            ax.set_title(f"Extreme values with {rp}-year return period \n based on {settings[f'yearmin1']}-{settings[f'yearmax1']}",fontsize=14)        

            ax = global_map(axs[1])
            bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[1,:]-rp_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); 
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=14)
            ax.set_title(f"Difference in extreme values with {rp}-year return period \n between {settings['yearmin1']}-{settings['yearmax1']} vs. {settings['yearmin0']}-{settings['yearmax0']} dataset",fontsize=14)        
=======
            bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:],transform=crt.crs.PlateCarree(),cmap=cmap3, vmin=vrange1[0], vmax=vrange1[1]); 
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=14)
            ax.set_title(f"Extreme values with {rp}-year return period \n based on {settings[f'yearmin{3}']}-{settings[f'yearmax{3}']}",fontsize=14)        

            ax = global_map(axs[1])
            bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:]-rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); 
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=14)
            ax.set_title(f"Difference in extreme values with {rp}-year return period \n between {settings['yearmin3']}-{settings['yearmax3']} vs. {settings['yearmin2']}-{settings['yearmax2']} dataset",fontsize=14)        
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            figname = f'EVA_map_RP{rp}_comparison_between_periods_v5.png' 
            fig.savefig(f'{dir_eva}/{figname}')

        # combined plot with 95th perc and 100 year return value
        for rp in [100]:
<<<<<<< HEAD
            rp_all = np.vstack((np.array(np.array(ds_gtsm_eva0[f'{str(rp)}_bf']),np.array(ds_gtsm_eva1[f'{str(rp)}_bf'])))
=======
            rp_all = np.vstack((np.array(ds_gtsm_eva1[f'{str(rp)}_bf']),np.array(ds_gtsm_eva2[f'{str(rp)}_bf']),np.array(ds_gtsm_eva3[f'{str(rp)}_bf'])))
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            
            mpl.rcParams.update({'font.size': 18})
            csize = 15
            fig = plt.figure(figsize=(20,15))
            axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
            axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 2), share_all=True, axes_pad=1.6,cbar_location='bottom',cbar_mode='each',cbar_size='8%',cbar_pad=0.4, label_mode='keep')
    
            # plot q95 and RP100 map
            ax = global_map(axs[0])
<<<<<<< HEAD
            bs = ax.scatter(x=coord_x,y=coord_y,s=csize,c=quan95[1,:],transform=crt.crs.PlateCarree(),cmap=cmap3, vmin=0, vmax=5);
=======
            bs = ax.scatter(x=coord_x,y=coord_y,s=csize,c=quan95[2,:],transform=crt.crs.PlateCarree(),cmap=cmap3, vmin=0, vmax=5);
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=20)
            ax.set_title('95th percentile values based on GTSM-ERA5-E',fontsize=22);

            ax = global_map(axs[1])
<<<<<<< HEAD
            bs = ax.scatter(x=coord_x,y=coord_y,s=csize,c=rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap3, vmin=0, vmax=8); 
=======
            bs = ax.scatter(x=coord_x,y=coord_y,s=csize,c=rp_all[2,:],transform=crt.crs.PlateCarree(),cmap=cmap3, vmin=0, vmax=8); 
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=20)
            ax.set_title('100-year return values based on GTSM-ERA5-E',fontsize=22);

            # plot q95 and RP100 bias
            ax = global_map(axs[2])
<<<<<<< HEAD
            bs =ax.scatter(x=coord_x,y=coord_y,s=csize,c=quan95[1,:]-quan95[0,:],cmap=cmap4,transform=crt.crs.PlateCarree(),vmin=-0.03, vmax=0.03);
=======
            bs =ax.scatter(x=coord_x,y=coord_y,s=csize,c=quan95[2,:]-quan95[1,:],cmap=cmap4,transform=crt.crs.PlateCarree(),vmin=-0.03, vmax=0.03);
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=20)
            ax.set_title('Difference in 95th percentile values \n GTSM-ERA5-E (1950-2022) vs. GTSM-ERA5 (1979-2018)',fontsize=22);

            ax = global_map(axs[3])
<<<<<<< HEAD
            bs =ax.scatter(x=coord_x,y=coord_y,s=csize,c=rp_all[1,:]-rp_all[0,:],cmap=cmap4,transform=crt.crs.PlateCarree(),vmin=-0.5, vmax=0.5); #
=======
            bs =ax.scatter(x=coord_x,y=coord_y,s=csize,c=rp_all[2,:]-rp_all[1,:],cmap=cmap4,transform=crt.crs.PlateCarree(),vmin=-0.5, vmax=0.5); #
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=20); 
            ax.set_title('Difference in 100-year return values \n GTSM-ERA5-E (1950-2022) vs. GTSM-ERA5 (1979-2018)',fontsize=22);
            plt.tight_layout()   

            figname = 'Map_comparison_model_GTSM-ERA5-E_vs_GTSM-ERA5.png'  
            fig.savefig(f'{dir_eva}/{figname}')


        # plot of difference in confidence interval width
        for rp in [100]:
<<<<<<< HEAD
            rp_bf_all = np.vstack((np.array(ds_gtsm_eva0[f'{str(rp)}_bf']),np.array(ds_gtsm_eva1[f'{str(rp)}_bf'])))
            rp_lower_all = np.vstack((np.array(ds_gtsm_eva0[f'{str(rp)}_lower']),np.array(ds_gtsm_eva1[f'{str(rp)}_lower'])))
            rp_higher_all = np.vstack((np.array(ds_gtsm_eva0[f'{str(rp)}_higher']),np.array(ds_gtsm_eva1[f'{str(rp)}_higher'])))
=======
            rp_bf_all = np.vstack((np.array(ds_gtsm_eva2[f'{str(rp)}_bf']),np.array(ds_gtsm_eva3[f'{str(rp)}_bf'])))
            rp_lower_all = np.vstack((np.array(ds_gtsm_eva2[f'{str(rp)}_lower']),np.array(ds_gtsm_eva3[f'{str(rp)}_lower'])))
            rp_higher_all = np.vstack((np.array(ds_gtsm_eva2[f'{str(rp)}_higher']),np.array(ds_gtsm_eva3[f'{str(rp)}_higher'])))
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            rp_intwidth_all = rp_higher_all - rp_lower_all; del rp_lower_all, rp_higher_all

            # count points where the width is reduced and the average reduction
            print(f'Confidence interval width reduced by {(1-np.mean(rp_intwidth_all[1,:]/rp_intwidth_all[0,:]))*100:.1f}%')
            print(f'Confidence interval width is reduced at {sum(1 for i in tmp if i < 0)/len(tmp)*100:.1f}% of points')
            
            mpl.rcParams.update({'font.size': 18})
            csize = 15
            fig = plt.figure(figsize=(20,15))
            axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
            axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(1, 2), share_all=True, axes_pad=1.6,cbar_location='bottom',cbar_mode='each',cbar_size='8%',cbar_pad=0.4, label_mode='keep')
    
            # plot interval width in the ERA5-E dataset
            ax = global_map(axs[0])
            bs = ax.scatter(x=coord_x,y=coord_y,s=csize,c=rp_intwidth_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap3, vmin=0, vmax=1);
            cbar = ax.cax.colorbar(bs); cbar.set_label('Confidence interval width [m]',fontsize=20)
            ax.set_title(f'Width of confidence interval \n {rp}-year return period water levels \n GTSM-ERA5-E dataset',fontsize=22);

            ax = global_map(axs[1])
            bs = ax.scatter(x=coord_x,y=coord_y,s=csize,c=rp_intwidth_all[1,:]-rp_intwidth_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=-0.5, vmax=0.5); 
            cbar = ax.cax.colorbar(bs); cbar.set_label('Difference [m]',fontsize=20)
            ax.set_title(f'Difference in confidence interval width \n {rp}-year return period water levels \n GTSM-ERA5-E (1950-2022) vs. GTSM-ERA5 (1979-2018)',fontsize=22);

            plt.tight_layout()   

            figname = 'Map_comparison_model_GTSM-ERA5-E_vs_GTSM-ERA5_confidence_intervals.png'  
            fig.savefig(f'{dir_eva}/{figname}')


    if make_local_plots:
    # Plotting comparison between periods with timeseries per location
    
        coords_locs = [[52.481, 4.506], # IJmuiden
                   [-4.0884, 39.734], # Mombasa, Kenya
                  [-34.381, -58.186], # Buenos Aires
                  [29.242, -94.767], # Houston
                  [40.504, -73.741], #NY
                  [31.574, 131.417], # Aburatsu, Japan
                  [42.341, -70.984], #Boston
                  [-22.923, -43.154], #Rio de J
                  [60.142, -1.145], #Lerwick, UK
                  [-6.163, 39.161]] # Zanzibar, Tanzania
        citynames = ['IJmuiden, NL','Mombasa, Kenya','Buenos Aires','Houston, USA','New York, USA','Aburatsu, Japan','Boston, USA','Rio de Janeiro, Brazil','Lerwick, UK','Zanzibar, Tanzania']
    
        dir_wlts = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly/waterlevel/'
        dir_wlts2 = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/'

        file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_1950_01_v1.nc')
        ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
        ds.load()

        loc_ind=[]
        for ii in range(0,10):
            abslat = np.abs(ds.station_y_coordinate-coords_locs[ii][0])
            abslon = np.abs(ds.station_x_coordinate-coords_locs[ii][1])
            loc_ind.append(np.argmin(abslon.values**2 + abslat.values**2))

        ids = loc_ind
        stations = ds.stations.values
    
        #load long timeseries
#        for year in range(1950,2023):
#            print(f'loading {year}')
#            for mnth in range(1,13):
#                print(f'loading {year}, month {mnth}')
#                if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
#                    file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                    ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#                    if "station_x_coordinate" in list(ds.data_vars):
#                    ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))                 
#                else:
#                    file_nc = os.path.join(dir_wlts2,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                    ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#                ds=ds.sel(stations=stations[ids],drop=True)
#                ds.load()
#                if ((year == 1950) & (mnth == 1)):
#                    ds_gtsm = ds
#                else:
#                    ds_gtsm = xr.concat([ds_gtsm,ds],dim="time")
#        # save timeseries selection 
#        ds_gtsm.to_netcdf("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1950_2022_1hr_selected_stations_CityLocs.nc");

        # open timeseries selection
        ds_gtsm = xr.open_dataset("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1950_2022_1hr_selected_stations_CityLocs.nc");
        
        ds_gtsm['sea_level'] = ds_gtsm['waterlevel']
        ds_gtsm = ds_gtsm.drop(['waterlevel'])
        ds_gtsm = detrend(ds_gtsm)
        ds_gtsm = ds_gtsm.drop(['sea_level'])
        ds_gtsm = ds_gtsm.chunk({"time": -1, "stations": "auto"})
        ds_gtsm.load()

        vrange1=[-0.5,0.5]; vrange2=[0,0.02]; 
        cmap = mpl.colormaps['seismic']#.resampled(20)

<<<<<<< HEAD
        ds_gtsm_eva0 = ds_gtsm_eva0.reset_index()
        ds_gtsm_eva1 = ds_gtsm_eva1.reset_index()
=======
        ds_gtsm_eva2 = ds_gtsm_eva2.reset_index()
        ds_gtsm_eva3 = ds_gtsm_eva3.reset_index()
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b

        overview_location_plot = 0
        comparison_location_plot = 1
    
    if overview_location_plot:
        # make overview plots per location
        for ss in range(0,len(ids)):
            
            st_id = ids[ss]
            cityname = citynames[ss]
            if stations[st_id] not in ds_gtsm['stations'].values:
                continue
            
            # prepare EVA outputs
            var = ds_gtsm.sea_level_detrended.sel(stations=stations[st_id]).to_dataframe().loc[:, 'sea_level_detrended'].dropna()
<<<<<<< HEAD
            varth = var.quantile(0.99); 
=======
            varth = var.quantile(0.995); 
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            model = EVA(var); del var
            model.get_extremes("POT", threshold=varth, r="72H");  del varth
            model.fit_model(distribution='genpareto',model='MLE')

            var = ds_gtsm.sea_level_detrended.sel(stations=stations[st_id],time=slice("1979-01-01", "2019-01-01")).to_dataframe().loc[:, 'sea_level_detrended'].dropna()
<<<<<<< HEAD
            varth = var.quantile(0.99)
=======
            varth = var.quantile(0.995)
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            model_short = EVA(var); del var
            model_short.get_extremes("POT", threshold=varth, r="72H"); del varth
            model_short.fit_model(distribution='genpareto',model='MLE')

            # get pandas with EVA results for plotting
            model_summary = model.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)
            model_short_summary = model_short.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)

            # plot
            fig = plt.figure(figsize=(20,20))
            ax0 = plt.subplot2grid((4, 2), (0, 0), colspan=1, projection=crt.crs.Robinson())
            ax1 = plt.subplot2grid((4, 2), (0, 1), colspan=1, rowspan=1)
            ax2 = plt.subplot2grid((4, 2), (1, 0), colspan=1, rowspan=1)
            ax3 = plt.subplot2grid((4, 2), (1, 1), colspan=1, rowspan=1,sharex=ax2,sharey=ax2)
            ax4 = plt.subplot2grid((4, 2), (2, 0), colspan=2, rowspan=1)
            ax5 = plt.subplot2grid((4, 2), (3, 0), colspan=2, rowspan=1)

        #    plot location on the map
            ax0 = global_map(ax0)
<<<<<<< HEAD
            bs0 = ax0.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[1,:]-rp_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
=======
            bs0 = ax0.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:]-rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
            bs = ax0.scatter(x=coord_x[st_id],y=coord_y[st_id],marker ='o',s=200,transform=crt.crs.PlateCarree(),facecolors='none',edgecolors='red',linewidth=3); 
            ax0.title.set_text(f"Location {int(stations[st_id])}")
            
        #   plot return values on a log plot per dataset
            model.plot_return_values(return_period=[1, 2, 5, 10, 25, 50, 100,500,1000],alpha=0.95,ax=ax2)
            model_short.plot_return_values(return_period=[1, 2, 5, 10, 25, 50, 100,500,1000],alpha=0.95,ax=ax3)
            ax2.set_ylabel('Still water level [m]'); ax3.set_ylabel('Still water level [m]')
            ax2.title.set_text(f"Return values based on 1950-2022 dataset"); ax3.title.set_text(f"Return values based on 1979-2018 dataset"); 

            # plot EVA confidence intervals together for comparison
            ax1.semilogx()
            ax1.grid(True, which="both")
            ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
            
            for col in ["lower ci", "upper ci"]:
                ax1.plot(model_summary.index.values,model_summary.loc[:, col].values,color="#5199FF",lw=1,ls="--",zorder=14)
                ax1.plot(model_short_summary.index.values,model_short_summary.loc[:, col].values,color="#ff6e51",lw=1,ls="--",zorder=15)
            ax1.fill_between(model_summary.index.values,model_summary.loc[:, "lower ci"].values,model_summary.loc[:, "upper ci"].values, facecolor="#5199FF",edgecolor="None",alpha=0.25,zorder=9)
            ax1.fill_between(model_short_summary.index.values,model_short_summary.loc[:, "lower ci"].values,model_short_summary.loc[:, "upper ci"].values, facecolor="#ff6e51",edgecolor="None",alpha=0.25,zorder=10)
            ax1.plot(model_summary.index.values,model_summary.loc[:, "return value"].values,color="#285fad",lw=2,ls="-",zorder=24,label='Based on 1950-2022 dataset')
            ax1.plot(model_short_summary.index.values,model_short_summary.loc[:, "return value"].values,color="#b8391f",lw=2,ls="-",zorder=25,label='Based on 1979-2018 dataset')
            ax1.set_xlabel("Return period")
            ax1.set_ylabel("Still water level [m]")
            ax1.legend(loc="upper left")
            ax1.title.set_text(f"Effect of timeseries length on modelled return values")

            # plot timeseries
            ts = ax4.plot(ds_gtsm.sel(stations=stations[st_id])['time'].values,ds_gtsm.sel(stations=stations[st_id])['sea_level_detrended'].values);
            ax4.set_ylabel('Still water level [m]')
            ax4.grid(); ax4.title.set_text('Full timeseries')

            # plot zoomed in timeseries
            peak_id = np.nanargmax(ds_gtsm.sel(stations=stations[st_id])['sea_level_detrended'].values)
            peak_time = ds_gtsm['time'].values[peak_id]
            time1 = peak_time - pd.DateOffset(days=3)
            time2 = peak_time + pd.DateOffset(days=3)
            ts = ax5.plot(ds_gtsm.sel(stations=stations[st_id]).sel(time=slice(time1,time2))['time'].values,ds_gtsm.sel(stations=stations[st_id]).sel(time=slice(time1,time2))['sea_level_detrended'].values,'bo--');
            ax5.set_ylabel('Still water level [m]')
            ax5.grid(); ax5.title.set_text('Timeseries around peak value')
            
            figname = f'EVA_station_{str(int(stations[st_id])).zfill(5)}_RV_1950-2022_vs_1979-2018_{cityname}.png' 
            fig.savefig(f'{dir_eva}/timeseries_plots/{figname}')
            plt.clf()

    if comparison_location_plot:
        # make comparison plot

        fig = plt.figure(figsize=(20,17))
        axs=[]
        axs.append(plt.subplot2grid((3, 3), (0, 0), colspan=1,rowspan=1))
        axs.append(plt.subplot2grid((3, 3), (0, 1), colspan=1, rowspan=1,sharex=axs[0]))
        axs.append(plt.subplot2grid((3, 3), (0, 2), colspan=1, rowspan=1,sharex=axs[0]))
        axs.append(plt.subplot2grid((3, 3), (1, 0), colspan=1, rowspan=1,sharex=axs[0]))
        axs.append(plt.subplot2grid((3, 3), (1, 1), colspan=1, rowspan=1,projection=crt.crs.Robinson()))
        axs.append(plt.subplot2grid((3, 3), (1, 2), colspan=1, rowspan=1,sharex=axs[0]))
        axs.append(plt.subplot2grid((3, 3), (2, 0), colspan=1, rowspan=1,sharex=axs[0]))
        axs.append(plt.subplot2grid((3, 3), (2, 1), colspan=1, rowspan=1,sharex=axs[0]))
        axs.append(plt.subplot2grid((3, 3), (2, 2), colspan=1, rowspan=1,sharex=axs[0]))

        cities = ['New York, USA','Lerwick, UK','IJmuiden, NL','Houston, USA','Aburatsu, Japan','Buenos Aires, Argentina','Rio de Janeiro, Brazil','Mombasa, Kenya']

        #add global map
<<<<<<< HEAD
        rp_diff = np.array(ds_gtsm_eva1['100_bf']) - np.array(ds_gtsm_eva0['100_bf'])
=======
        rp_diff = np.array(ds_gtsm_eva3['100_bf']) - np.array(ds_gtsm_eva2['100_bf'])
>>>>>>> 2682cc3d03a9e92d2a7f309474ed169b67916d9b
        axs[4] = global_map(axs[4])
        bs0 = axs[4].scatter(x=coord_x,y=coord_y,s=15,c=rp_diff,transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
        axs[4].title.set_text(f"Difference in water level \n 100-year return values")
        cbar = plt.colorbar(bs0,ax=axs[4],orientation='horizontal',aspect=30); cbar.set_label('Still water level difference [m]')

        for cc in range(0,len(cities)):
            cname = cities[cc]

            if cc<4:
                ax = axs[cc]
            else:
                ax = axs[cc+1]

            st_id = ids[citynames.index(cname)]

            if stations[st_id] not in ds_gtsm['stations'].values:
                continue

            # plot location on the map
            bs = axs[4].scatter(x=coord_x[st_id],y=coord_y[st_id],marker ='o',s=150,transform=crt.crs.PlateCarree(),facecolors='none',edgecolors='red',linewidth=3); 
            
            # prepare EVA outputs
            var = ds_gtsm.sea_level_detrended.sel(stations=stations[st_id]).to_dataframe().loc[:, 'sea_level_detrended'].dropna()
            varth = var.quantile(0.99); 
            model = EVA(var); del var
            model.get_extremes("POT", threshold=varth, r="72H");  del varth
            model.fit_model(distribution='genpareto',model='MLE')

            var = ds_gtsm.sea_level_detrended.sel(stations=stations[st_id],time=slice("1979-01-01", "2019-01-01")).to_dataframe().loc[:, 'sea_level_detrended'].dropna()
            varth = var.quantile(0.99)
            model_short = EVA(var); del var
            model_short.get_extremes("POT", threshold=varth, r="72H"); del varth
            model_short.fit_model(distribution='genpareto',model='MLE')

            # get pandas with EVA results for plotting
            model_summary = model.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)
            model_short_summary = model_short.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)

            # plot EVA confidence intervals together for comparison
            ax.semilogx()
            ax.grid(True, which="both",color='darkgrey')
            ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
            
            for col in ["lower ci", "upper ci"]:
                ax.plot(model_summary.index.values,model_summary.loc[:, col].values,color="#5199FF",lw=1,ls="--",zorder=14)
                ax.plot(model_short_summary.index.values,model_short_summary.loc[:, col].values,color="#ff6e51",lw=1,ls="--",zorder=15)
            ax.fill_between(model_summary.index.values,model_summary.loc[:, "lower ci"].values,model_summary.loc[:, "upper ci"].values, facecolor="#5199FF",edgecolor="None",alpha=0.25,zorder=9)
            ax.fill_between(model_short_summary.index.values,model_short_summary.loc[:, "lower ci"].values,model_short_summary.loc[:, "upper ci"].values, facecolor="#ff6e51",edgecolor="None",alpha=0.25,zorder=10)
            ax.plot(model_summary.index.values,model_summary.loc[:, "return value"].values,color="#285fad",lw=2,ls="-",zorder=24,label='GTSM-ERA5-E (1950-2022)')
            ax.plot(model_short_summary.index.values,model_short_summary.loc[:, "return value"].values,color="#b8391f",lw=2,ls="-",zorder=25,label='GTSM-ERA5 (1979-2018)')

            if cc>4:
                ax.set_xlabel("Return period")
            if cc in [0,3,5]:
                ax.set_ylabel("Still water level [m]")
            if cc==0:
                ax.legend(loc="upper left")
            ax.title.set_text(cname)
            plt.tight_layout()
            
            figname = f'EVA_station_CITIES_RV_1950-2022_vs_1979-2018.png' 
            fig.savefig(f'{dir_eva}/{figname}')



        

