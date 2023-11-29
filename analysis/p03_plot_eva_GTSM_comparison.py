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
# Contact: natalia.aleksandrova@deltares.nl
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
import matplotlib.ticker as tkr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import cartopy as crt
from eva import detrend
from pyextremes import EVA
from pyextremes.plotting import plot_return_values

#from pyextremes import plot_threshold_stability

if __name__ == "__main__":   
    settings = {'yearmin1': 1985,'yearmax1': 2014,'mode1':'1hr',
                'yearmin2': 1979,'yearmax2': 2018,'mode2':'1hr',
                'yearmin3': 1950,'yearmax3': 2022,'mode3':'1hr'}
    rps =[1,10,50,100]    
    
    # location of EVA data
    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_eva = os.path.join(dir_postproc,'EVA-GTSM-ERA5')
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
    file_list_nc = [str(file_list1[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list1))]
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
    file_list_nc = [str(file_list1[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list1))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    quan90_1 = ds['sea_level_detrended'].values[5,:]
    quan95_1 = ds['sea_level_detrended'].values[6,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        quan90_1 = np.append(quan90_1,ds['sea_level_detrended'].values[5,:])
        quan95_1 = np.append(quan95_1,ds['sea_level_detrended'].values[6,:])
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
    del ds, file_list_nc    

    # loading CDS return values and percentiles
#    dir_rv = os.path.join(dir_eva,'EVA_RV_from_CDS')
#    cds_rv100 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp100_best-fit_v1.nc'); cds_rv100.close()
#    cds_rv50 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp50_best-fit_v1.nc'); cds_rv50.close()
#    cds_rv10 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp10_best-fit_v1.nc'); cds_rv10.close()
#    cds_rv1 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp1_best-fit_v1.nc'); cds_rv1.close()

#    cds_rv = np.array([cds_rv1.sel(stations=ds_gtsm_eva1['station'].values)['return_mean_water_level'].values,
#                     cds_rv10.sel(stations=ds_gtsm_eva1['station'].values)['return_mean_water_level'].values,
#                     cds_rv50.sel(stations=ds_gtsm_eva1['station'].values)['return_mean_water_level'].values,
#                     cds_rv100.sel(stations=ds_gtsm_eva1['station'].values)['return_mean_water_level'].values])

#    cds_perc95 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_95-percentile_v1.nc'); cds_perc95.close()
#    cds_perc90 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_90-percentile_v1.nc'); cds_perc90.close()    

    # PLOTTING
    print('Plotting... ')
    cmap = mpl.colormaps['viridis'].resampled(20)
    cmap2 = mpl.colormaps['hot_r']#.resampled(20)
    cmap3 = mpl.colormaps['magma_r'].resampled(20)
    cmap4 = mpl.colormaps['seismic']#.resampled(20)
    
    # # Plotting percentiles - periods 1,2,3 and difference
    # quan95 = np.vstack((quan95_1,quan95_2,quan95_3))
    # quan90 = np.vstack((quan90_1,quan90_2,quan90_3))
    # vrange1=[0,3]; vrange2=[-0.03,0.03];   
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

    # # ---- Smaller plot with just two panes
    # fig = plt.figure(figsize=(20,10))
    # axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    # axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(1, 2), share_all=True, axes_pad=1.7,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=0.3, label_mode='keep')
    # ax = global_map(axs[0])
    # bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[2,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
    # cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=12)
    # ax.title.set_text(f"95th percentile of still water levels \n based on {settings[f'yearmin3']}-{settings[f'yearmax3']} timeseries")
    # ax = global_map(axs[1])
    # bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan95[2,:]-quan95[1,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); 
    # cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=12)
    # ax.title.set_text(f"Difference in 95th percentile values \n between {settings['yearmin3']}-{settings['yearmax3']} vs. {settings['yearmin2']}-{settings['yearmax2']} timeseries")
    # figname = 'EVA_map_95percentile_comparison_between_periods_v5.png'  
    # fig.savefig(f'{dir_eva}/{figname}')


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

    # # smaller plot with just two panels
    # for rp in [1,10,100]:
    #     rp_all = np.vstack((np.array(ds_gtsm_eva1[f'{str(rp)}_bf']),np.array(ds_gtsm_eva2[f'{str(rp)}_bf']),np.array(ds_gtsm_eva3[f'{str(rp)}_bf'])))
    #     vrange1=[0,3]; vrange2=[-0.5,0.5];   
    #     fig = plt.figure(figsize=(20,10))
    #     axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    #     axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(1, 2), share_all=True, axes_pad=1.7,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=0.3, label_mode='keep')

    #     ax = global_map(axs[0])
    #     bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
    #     cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=12)
    #     ax.title.set_text(f"Extreme value for {rp}-year return period \n based on {settings[f'yearmin{3}']}-{settings[f'yearmax{3}']}")        

    #     ax = global_map(axs[1])
    #     bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:]-rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap4, vmin=vrange2[0], vmax=vrange2[1]); 
    #     cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level difference [m]',fontsize=12)
    #     ax.title.set_text(f"Difference in extreme values for {rp}-year return period \n between {settings['yearmin3']}-{settings['yearmax3']} vs. {settings['yearmin2']}-{settings['yearmax2']} timeseries")        
    #     figname = f'EVA_map_RP{rp}_comparison_between_periods_v5.png' 
    #     fig.savefig(f'{dir_eva}/{figname}')

    



    # Plotting comparison between periods with timeseries per location
    
    rp_all = np.vstack((np.array(ds_gtsm_eva1['100_bf']),np.array(ds_gtsm_eva2['100_bf']),np.array(ds_gtsm_eva3['100_bf'])))
    diff100 = abs(rp_all[1,:]-rp_all[2,:])
    stations = np.array(ds_gtsm_eva1['station'])
    ids = [i for i,v in enumerate(diff100) if v > 0.5]
    ids0= [i for i in np.arange(0,43000,1000)]
    ids = ids +ids0

    diff100_sel = diff100[ids]; coordx_sel = coord_x[ids]; coordy_sel = coord_y[ids]; st_sel = stations[ids]

    dir_wlts = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly/waterlevel/'
    dir_wlts2 = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/'
    
    #load timeseries
#    for year in range(1950,2022):
#        print(f'loading {year}')
#        for mnth in range(1,13):
#            print(f'loading {year}, month {mnth}')
#            if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
#                file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#                if "station_x_coordinate" in list(ds.data_vars):
#                    ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))                 
#            else:
#                file_nc = os.path.join(dir_wlts2,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#            ds=ds.sel(stations=stations[ids],drop=True)
#            ds.load()
#            if ((year == 1950) & (mnth == 1)):
#                ds_gtsm = ds
#            else:
#                ds_gtsm = xr.concat([ds_gtsm,ds],dim="time")
#    ds_gtsm.to_netcdf("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1950_2022_1hr_selected_stations.nc");
    ds_gtsm = xr.open_dataset("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1950_2022_1hr_selected_stations.nc");
    ds_gtsm['sea_level'] = ds_gtsm['waterlevel']
    ds_gtsm = ds_gtsm.drop(['waterlevel'])
    ds_gtsm = detrend(ds_gtsm)
    ds_gtsm = ds_gtsm.drop(['sea_level'])
    ds_gtsm = ds_gtsm.chunk({"time": -1, "stations": "auto"})
    ds_gtsm.load()

#    for year in range(1979,2023):
#        for mnth in range(1,13):
#            if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
#                file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#                if "station_x_coordinate" in list(ds.data_vars):
#                    ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))                 
#            else:
#                file_nc = os.path.join(dir_wlts2,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#            if ((year == 1950) & (mnth == 1)):
#                ds_gtsm_short = ds.sel(stations=stations[ids])
#            else:
#                ds_gtsm_short = xr.concat([ds_gtsm,ds.sel(stations=stations[ids])],dim="time")
#    ds_gtsm_short.to_netcdf("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1979_2022_1hr_selected_stations.nc");
#    ds_gtsm_short = xr.open_dataset("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1979_2022_1hr_selected_stations.nc");
#    ds_gtsm_short['sea_level'] = ds_gtsm_short['waterlevel']
#    ds_gtsm_short = ds_gtsm_short.drop(['waterlevel'])
#    ds_gtsm_short = detrend(ds_gtsm_short)
#    ds_gtsm_short = ds_gtsm_short.drop(['sea_level'])
#    ds_gtsm_short = ds_gtsm_short.chunk({"time": -1, "stations": "auto"})
#    ds_gtsm_short.load()

    vrange1=[-0.5,0.5]; vrange2=[0,0.02]; 
    cmap = mpl.colormaps['seismic']#.resampled(20)

    ds_gtsm_eva2 = ds_gtsm_eva2.reset_index()
    ds_gtsm_eva3 = ds_gtsm_eva3.reset_index()
    
    # make overview plot
    for st_id in ids:
        #station = stations[st_id] 

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
        model_summary = model.get_summary(return_period=[1,2,5,10,50,100],alpha=0.95)
        model_short_summary = model_short.get_summary(return_period=[1,2,5,10,50,100],alpha=0.95)

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
        bs0 = ax0.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[2,:]-rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
        bs = ax0.scatter(x=coord_x[st_id],y=coord_y[st_id],marker ='o',s=200,transform=crt.crs.PlateCarree(),facecolors='none',edgecolors='red',linewidth=3); 
        ax0.title.set_text(f"Location {int(stations[st_id])}")
        
    #   plot return values on a log plot per dataset
        model.plot_return_values(return_period=[1, 2, 5, 10, 25, 50, 100],alpha=0.95,ax=ax2)
        model_short.plot_return_values(return_period=[1, 2, 5, 10, 25, 50, 100],alpha=0.95,ax=ax3)
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
        ax1.title.set_text(f"Effect of timeseries length on return values")

        # plot timeseries
        ts = ax4.plot(ds_gtsm.sel(stations=stations[st_id])['time'].values,ds_gtsm.sel(stations=stations[st_id])['sea_level_detrended'].values);
        ax4.set_ylabel('Still water level [m]')
        ax4.grid(); ax4.title.set_text('Full timeseries')

        # plot zoomed in timeseries
        peak_id = np.nanargmax(ds_gtsm.sel(stations=stations[st_id])['sea_level_detrended'].values)
        peak_time = ds_gtsm['time'].values[peak_id]
        time1 = peak_time - pd.DateOffset(days=3)
        time2 = peak_time + pd.DateOffset(days=3)
        ts = ax5.plot(ds_gtsm.sel(stations=stations[st_id]).sel(time=slice(time1,time2))['time'].values,ds_gtsm.sel(stations=stations[st_id]).sel(time=slice(time1,time2))['sea_level_detrended'].values);
        ax5.set_ylabel('Still water level [m]')
        ax5.grid(); ax5.title.set_text('Timeseries around peak value')
        
        figname = f'EVA_station_{str(int(stations[st_id])).zfill(5)}_RV_1950-2022_vs_1979-2018.png' 
        fig.savefig(f'{dir_eva}/timeseries_plots/{figname}')
        plt.clf()


    

    # COMPARISON WITH CDS VALUES
    
    #check RVs for model and CDS with timeseries
#    vrange1=[0,3]; vrange2=[0,0.02]; 
#    cmap = mpl.colormaps['viridis']
    #dir_wlts = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly/waterlevel/'
#    dir_wlts2 = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/'

#    locs = np.arange(0,44)*1000+50
    
    #load timeseries
#    for year in range(1985,2014):
#        for mnth in range(1,13):
#            if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
#                file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#                if "station_x_coordinate" in list(ds.data_vars):
#                    ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate")) 
#                if "year" in list(ds.coords):
#                    ds = ds.drop("year")
#            else:
#                file_nc = os.path.join(dir_wlts2,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#            if ((year == 1985) & (mnth == 1)):
#                ds_gtsm = ds.sel(stations=locs)
#            else:
#                ds_gtsm = xr.concat([ds_gtsm,ds.sel(stations=locs)],dim="time")

    # detrend
#    ds_gtsm['sea_level'] = ds_gtsm['waterlevel']
#    ds_gtsm = ds_gtsm.drop(['waterlevel'])
#    ds_gtsm = detrend(ds_gtsm)
#    ds_gtsm = ds_gtsm.drop(['sea_level'])
#    ds_gtsm = ds_gtsm.chunk({"time": -1, "stations": "auto"})
#    ds_gtsm.load()
#
#    for st in locs:
#        st_id = np.where(stations==st); st_id = st_id[0][0]
        
#        fig = plt.figure(figsize=(20,10))
#        ax0 = plt.subplot2grid((2, 2), (0, 0), colspan=1, projection=crt.crs.Robinson())
#        ax1 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)
#        ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2, rowspan=1)

    #    plot location on the map
#        ax0 = global_map(ax0)
#        bs0 = ax0.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
#        bs = ax0.scatter(x=coord_x[st_id],y=coord_y[st_id],s=80,c='red',transform=crt.crs.PlateCarree(),facecolors='none', edgecolors='r'); 
#        ax0.title.set_text(f"Location {stations[st_id]}")
        
    #   plot return values on a log plot
#        eva_data = ds_gtsm_eva1.loc[(ds_gtsm_eva1['station'].values==stations[st_id])]
#        tmp1 = ax1.scatter([1,10,50,100],eva_data[['1_bf','10_bf','50_bf','100_bf']],label=f"Current model data, {settings['yearmin1']}-{settings['yearmax1']}"); ax1.semilogx(); plt.grid();
#        tmp2 = ax1.scatter([1,10,50,100],cds_rv[:,st_id],c='red',label=f"CDS, {settings['yearmin2']}-{settings['yearmax2']}");
#        ax1.legend(loc="upper left")
        #ax1.set_xlabel('Return period in years')
#        ax1.set_ylabel('Total water level [m]')
#        ax1.grid()
#        ax1.title.set_text(f"Return values"); 
        
        # plot timeseries
#        ts = ax2.plot(ds_gtsm.sel(stations=stations[st_id])['time'].values,ds_gtsm.sel(stations=stations[st_id])['sea_level_detrended'].values);
#        ax2.set_ylabel('Total water level [m]')
#        
#        figname = f'EVA_station_{str(stations[st]).zfill(5)}_RV_CDSvsModel_with_timeseries.png' 
#        fig.savefig(f'{dir_eva}/{figname}')
#        plt.clf(fig)
        
    # # Plot alternative using pyextremes
    
    # for st_id in ids[0::50]: 

    #     # POT data
    #     q=0.998
    #     rps = [1, 2, 5, 10, 25, 50, 100]
    #     var1 = ds_gtsm.sea_level_detrended.sel(stations=stations[st_id])
    #     var1 = var1.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
    #     varth = var1.quantile(q)
    #     model1 = EVA(var1)
    #     model1.get_extremes("POT", threshold=varth, r="72H"); #model.plot_extremes(show_clusters=True)
    #     model1.fit_model()  
    #     summary1 = model1.get_summary(return_period=rps,alpha=0.95)
        
    #     var2 = ds_gtsm.sel(time=slice('1979-01-01','2023-01-01'))
    #     var2 = var2.sea_level_detrended.sel(stations=stations[st_id])
    #     var2 = var2.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
    #     varth = var2.quantile(q)
    #     model2 = EVA(var2)
    #     model2.get_extremes("POT", threshold=varth, r="72H")
    #     model2.plot_extremes(show_clusters=True)
    #     model2.fit_model()
    #     summary2 = model2.get_summary(return_period=rps,alpha=0.95)

    #     fig = plt.figure(figsize=(20,10))
    #     ax0 = plt.subplot2grid((2, 2), (0, 0), colspan=1, projection=crt.crs.Robinson())
    #     ax1 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)
    #     ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=1)
    #     ax3 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=1)

    #     # plot location on the map
    #     ax0 = global_map(ax0)
    #     bs0 = ax0.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
    #     bs = ax0.scatter(x=coord_x[st_id],y=coord_y[st_id],s=80,c='red',transform=crt.crs.PlateCarree(),facecolors='none', edgecolors='r'); 
    #     ax0.title.set_text(f"Location {stations[st_id]}")
        
    #     # plot return values on a log plot
    #     tmp1 = ax1.scatter(rps,summary1['return value'],label=f"{settings['yearmin1']}-{settings['yearmax1']}"); ax1.semilogx(); plt.grid();
    #     tmp2 = ax1.scatter(rps,summary2['return value'],c='red',label=f"{settings['yearmin2']}-{settings['yearmax2']}");
    #     ax1.legend(loc="upper left")
    #     #ax1.set_xlabel('Return period in years')
    #     ax1.set_ylabel('Total water level [m]')
    #     ax1.grid()
    #     ax1.title.set_text(f"Return values"); 
        
    #     # plot pot plot full TS
    #     model1.plot_return_values(return_period=np.logspace(0.01, 2, 100), alpha=0.95,ax=ax2); 
    #     ax2.title.set_text(f"Return values, {settings['yearmin1']}-{settings['yearmax1']}")
    #     plt.tight_layout()
        
    #     # plot pot plot part TS 
    #     model2.plot_return_values(return_period=np.logspace(0.01, 2, 100), alpha=0.95,ax=ax3); 
    #     ax3.title.set_text(f"Return values, {settings['yearmin2']}-{settings['yearmax2']}")
    #     plt.tight_layout()
    #     plt.show()
    #     figname = f'EVA_station_{str(stations[st_id]).zfill(5)}_pot_plots_pyextremes2.png' 
    #     fig.savefig(f'{dir_eva}/checkplots/{figname}')    

    #Plotting difference of present calculation and dataset on CDS - rv values - difference
#    vrange=[-0.5,0.5]    
#    print('Plotting... ')
#    fig = plt.figure(figsize=(26,20))
#    axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
#    axs = AxesGrid(fig, 111, axes_class=axes_class,
#           nrows_ncols=(2, 2),
#           share_all=True,
#           axes_pad=1.7,
#           cbar_location='right',
#           cbar_mode='each',
#           cbar_size='3%',
#           cbar_pad=1.0,
#           label_mode='keep')
#    for ii in range(0,4):
#       ax = global_map(axs[ii])
#       bs = ax.scatter(x=coord_x,y=coord_y,
#                       s=25,c=ds_gtsm_eva1[f'{str(rps[ii])}_bf'].values-cds_rv[ii,:],transform=crt.crs.PlateCarree(),
#                       cmap=cmap4, vmin=vrange[0], vmax=vrange[1])
#       ax.title.set_text(f"Total water level RP{rps[ii]}, difference model vs. CDS values \n period: {settings['yearmin1']}-{settings['yearmax1']}")
#       cbar = ax.cax.colorbar(bs)
#    figname = 'EVA_map_compare_dataset_to_CDS_1985-2014_RVdiff_v2.png' 
#    fig.savefig(f'{dir_eva}/{figname}')

    #Plotting difference of present calculation and dataset on CDS - rv values
 #   vrange=[0,3]    
 #   print('Plotting... ')
 #   fig = plt.figure(figsize=(26,20))
 #   axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
 #   axs = AxesGrid(fig, 111, axes_class=axes_class,
 #          nrows_ncols=(2, 2),
 #          share_all=True,
 #          axes_pad=1.7,
 #          cbar_location='right',
 #          cbar_mode='each',
 #          cbar_size='3%',
 #          cbar_pad=1.0,
 #          label_mode='keep')
    
 #   ax = global_map(axs[0])
 #   bs = ax.scatter(x=coord_x,y=coord_y,
 #                      s=25,c=ds_gtsm_eva1['10_bf'].values,transform=crt.crs.PlateCarree(),
 #                      cmap=cmap, vmin=vrange[0], vmax=vrange[1])
 #   ax.title.set_text(f"Total water level RP1, Model result for period: {settings['yearmin1']}-{settings['yearmax1']}"); cbar = ax.cax.colorbar(bs)

#    ax = global_map(axs[1])
#    bs = ax.scatter(x=coord_x,y=coord_y,
#                       s=25,c=cds_rv[1,:],transform=crt.crs.PlateCarree(),
#                       cmap=cmap, vmin=vrange[0], vmax=vrange[1])
#    ax.title.set_text(f"Total water level RP1, CDS data"); cbar = ax.cax.colorbar(bs)

#    ax = global_map(axs[2])
#    bs = ax.scatter(x=coord_x,y=coord_y,
#                    s=25,c=ds_gtsm_eva1['100_bf'].values - ds_gtsm_eva1['10_bf'].values,transform=crt.crs.PlateCarree(),
#                    cmap=cmap2, vmin=0, vmax=1)
#    ax.title.set_text(f"Difference RP100 vs. RP10 , Model result for period: {settings['yearmin1']}-{settings['yearmax1']}"); cbar = ax.cax.colorbar(bs)
#    #bs = ax.scatter(x=coord_x,y=coord_y,
    #                   s=25,c=ds_gtsm_eva1['100_bf'].values,transform=crt.crs.PlateCarree(),
    #                   cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    #ax.title.set_text(f"Total water level RP100, Model result for period: {settings['yearmin1']}-{settings['yearmax1']}"); cbar = ax.cax.colorbar(bs)

#    ax = global_map(axs[3])
#    bs = ax.scatter(x=coord_x,y=coord_y,
 #                      s=25,c=cds_rv[3,:]-cds_rv[1,:],transform=crt.crs.PlateCarree(),
  #                     cmap=cmap2, vmin=0, vmax=1)
    #ax.title.set_text(f"Difference RP100 vs. RP10, CDS data"); cbar = ax.cax.colorbar(bs)

#    figname = 'EVA_map_compare_dataset_to_CDS_1985-2014_RV_relative2_v2.png' 
#    fig.savefig(f'{dir_eva}/{figname}')    

