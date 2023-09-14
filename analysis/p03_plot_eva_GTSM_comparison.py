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
    settings = {'yearmin1': 1951,'yearmax1': 2022,'mode1':'1hr',
                'yearmin2': 1979,'yearmax2': 2022,'mode2':'1hr',
                'yearmin3': 1985,'yearmax3': 2014,'mode3':'10min'}
    rps =[2,10,50,100]    
    
    # location of EVA data
    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']
    dir_eva = os.path.join(dir_postproc,'EVA-GTSM-ERA5')
    dir_eva1 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin1"]}_{settings["yearmax1"]}_{settings["mode1"]}') #output dir
    dir_eva2 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin2"]}_{settings["yearmax2"]}_{settings["mode2"]}')
    dir_eva3 = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin3"]}_{settings["yearmax3"]}_{settings["mode3"]}')
    
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
    
    #locate and read files that contain coordinates 
    file_list_nc = [str(file_list1[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list1))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    coord_x = ds['station_x_coordinate'].values
    coord_y = ds['station_y_coordinate'].values
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        coord_x = np.append(coord_x, ds['station_x_coordinate'].values)
        coord_y = np.append(coord_y, ds['station_y_coordinate'].values)

    # read quantiles
    file_list_nc = [str(file_list1[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list1))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    #quan90_1 = ds['sea_level_detrended'].values[0,:]
    quan95_1 = ds['sea_level_detrended'].values[1,:]
    quan999_1 = ds['sea_level_detrended'].values[3,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        #quan90_1 = np.append(quan90_1,ds['sea_level_detrended'].values[0,:])
        quan95_1 = np.append(quan95_1,ds['sea_level_detrended'].values[1,:])
        quan999_1 = np.append(quan999_1,ds['sea_level_detrended'].values[3,:])
    del ds, file_list_nc
    
    file_list_nc = [str(file_list2[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list2))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    #quan90_2 = ds['sea_level_detrended'].values[0,:]
    quan95_2 = ds['sea_level_detrended'].values[1,:]
    quan999_2 = ds['sea_level_detrended'].values[3,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        #quan90_2 = np.append(quan90_2,ds['sea_level_detrended'].values[0,:])
        quan95_2 = np.append(quan95_2,ds['sea_level_detrended'].values[1,:])
        quan999_2 = np.append(quan999_2,ds['sea_level_detrended'].values[3,:])
    del ds, file_list_nc

    file_list_nc = [str(file_list3[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list3))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    #quan90_3 = ds['sea_level_detrended'].values[0,:]
    quan95_3 = ds['sea_level_detrended'].values[1,:]
    quan999_3 = ds['sea_level_detrended'].values[3,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        #quan90_3 = np.append(quan90_3,ds['sea_level_detrended'].values[0,:])
        quan95_3 = np.append(quan95_3,ds['sea_level_detrended'].values[1,:])
        quan999_3 = np.append(quan999_3,ds['sea_level_detrended'].values[3,:])
    del ds, file_list_nc    

    # loading CDS return values and percentiles
    dir_rv = os.path.join(dir_eva,'EVA_RV_from_CDS')
    cds_rv100 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp100_best-fit_v1.nc'); cds_rv100.close()
    cds_rv50 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp50_best-fit_v1.nc'); cds_rv50.close()
    cds_rv10 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp10_best-fit_v1.nc'); cds_rv10.close()
    cds_rv2 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_rp2_best-fit_v1.nc'); cds_rv2.close()

    cds_rv = np.array([cds_rv2.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values,
                    cds_rv10.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values,
                     cds_rv50.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values,
                    cds_rv100.sel(stations=ds_gtsm_eva3['station'].values)['return_mean_water_level'].values])

    cds_perc95 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_95-percentile_v1.nc'); cds_perc95.close()
    cds_perc90 = xr.open_dataset(f'{dir_rv}/reanalysis_waterlevel_actual-value_1985-2014_90-percentile_v1.nc'); cds_perc90.close()    

    # PLOTTING
    print('Plotting... ')
    cmap = mpl.colormaps['viridis'].resampled(20)
    cmap2 = mpl.colormaps['hot_r'].resampled(20)

    # Plotting percentiles - periods 1,2,3 and difference
    #quan95 = np.vstack((quan95_1,quan95_2,quan95_3))
    quan999 = np.vstack((quan999_1,quan999_2,quan999_3))
    vrange1=[0,3]; vrange2=[0,0.02];   
    fig = plt.figure(figsize=(26,20))
    axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 2), share_all=True, axes_pad=1.7,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=1.0, label_mode='keep')
    for ii in range(0,3):
        ax = global_map(axs[ii])
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan999[ii,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); cbar = ax.cax.colorbar(bs); 
        ax.title.set_text(f"99.9th perc., period {settings[f'yearmin{ii+1}']}-{settings[f'yearmax{ii+1}']}")
    ax = global_map(axs[3])
    bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=quan999[0,:]-quan999[2,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange2[0], vmax=vrange2[1]); cbar = ax.cax.colorbar(bs); 
    ax.title.set_text(f"99.9th perc. difference, {settings['yearmin1']}-{settings['yearmax1']} vs. {settings['yearmin3']}-{settings['yearmax3']}")
    figname = 'EVA_map_999percentile_comparison_between_periods.png' 
    fig.savefig(f'{dir_eva}/{figname}')

    # Plotting return values
    for rp in [5,10,100]:
        rp_all = np.vstack((np.array(ds_gtsm_eva1[str(rp)]),np.array(ds_gtsm_eva2[str(rp)]),np.array(ds_gtsm_eva3[str(rp)])))
        vrange1=[0,3]; vrange2=[0,0.5];   
        fig = plt.figure(figsize=(26,20))
        axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
        axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 2), share_all=True, axes_pad=1.7,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=1.0, label_mode='keep')
        ax = global_map(axs[0])
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); cbar = ax.cax.colorbar(bs); 
        ax.title.set_text(f"Total water level RP{rp}, period {settings[f'yearmin{1}']}-{settings[f'yearmax{1}']}")
        ax = global_map(axs[1])
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[1,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); cbar = ax.cax.colorbar(bs); 
        ax.title.set_text(f"Total water level RP{rp}, period {settings[f'yearmin{2}']}-{settings[f'yearmax{2}']}")        
        #bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=abs(rp_all[0,:]-rp_all[2,:]),transform=crt.crs.PlateCarree(),cmap=cmap2, vmin=vrange2[0], vmax=vrange2[1]); cbar = ax.cax.colorbar(bs); 
        #ax.title.set_text(f"Total water level RP{rp}, diff. {settings['yearmin1']}-{settings['yearmax1']} vs. {settings['yearmin3']}-{settings['yearmax3']}")
        ax = global_map(axs[2])
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=abs(rp_all[0,:]-rp_all[1,:]),transform=crt.crs.PlateCarree(),cmap=cmap2, vmin=vrange2[0], vmax=vrange2[1]); cbar = ax.cax.colorbar(bs); 
        ax.title.set_text(f"Total water level RP{rp}, diff. {settings['yearmin1']}-{settings['yearmax1']} vs. {settings['yearmin2']}-{settings['yearmax2']}")
        ax = global_map(axs[3])
        bs = ax.scatter(x=coord_x,y=coord_y,s=15,c=abs(rp_all[0,:]-rp_all[2,:]),transform=crt.crs.PlateCarree(),cmap=cmap2, vmin=vrange2[0], vmax=vrange2[1]); cbar = ax.cax.colorbar(bs); 
        ax.title.set_text(f"Total water level RP{rp}, diff. {settings['yearmin1']}-{settings['yearmax1']} vs. {settings['yearmin3']}-{settings['yearmax3']}")        
        figname = f'EVA_map_RP{rp}_comparison_between_periods.png' 
        fig.savefig(f'{dir_eva}/{figname}')

    # locate points in deep water where large differences are observed
    rp_all = np.vstack((np.array(ds_gtsm_eva1['100']),np.array(ds_gtsm_eva2['100']),np.array(ds_gtsm_eva3['100'])))
    diff100 = abs(rp_all[0,:]-rp_all[1,:])
    stations = np.array(ds_gtsm_eva1['station'])
    ids = [i for i,v in enumerate(diff100) if v > 0.5]
    diff100_sel = diff100[ids]; coordx_sel = coord_x[ids]; coordy_sel = coord_y[ids]; st_sel = stations[ids]
    
    # make a plot for a station in the Pacific with a large difference
    #ids2 = np.where((coordy_sel<12) & (coordy_sel>0) & (coordx_sel<-110) & (coordx_sel>-140))
    #stations1 = st_sel[ids2] # station numbers of the points of interest
    #id_st1 = np.squeeze(np.where(np.isin(stations,stations1)))
    #rp_all[0,id_st1]
    #rp_all[1,id_st1]
    #coord_x[id_st1]
    #coord_y[id_st1]

    vrange1=[0,3]; vrange2=[0,0.02]; 
    cmap = mpl.colormaps['viridis'].resampled(20)
    dir_wlts = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly/waterlevel/'
    dir_wlts2 = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/'
    
    #load timeseries
#    for year in range(1951,2023):
#        for mnth in range(1,13):
#            if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
#               ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#                if "station_x_coordinate" in list(ds.data_vars):
#                    ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))                 
#            else:
#                file_nc = os.path.join(dir_wlts2,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
#                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
#            if ((year == 1951) & (mnth == 1)):
#                ds_gtsm = ds.sel(stations=stations[ids])
#            else:
#                ds_gtsm = xr.concat([ds_gtsm,ds.sel(stations=stations[ids])],dim="time")
#    ds_gtsm.to_netcdf("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1951_2022_1hr_selection_timeseries.nc");
ds_gtsm = xr.open_dataset("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1951_2022_1hr_selection_timeseries.nc");

for st_id in ids[0::100]: 
        fig = plt.figure(figsize=(26,5))
        ax0 = plt.subplot2grid((2, 2), (0, 0), colspan=1, projection=crt.crs.Robinson())
        ax1 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)
        ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2, rowspan=1)
        ax0 = global_map(ax0)
        bs0 = ax0.scatter(x=coord_x,y=coord_y,s=15,c=rp_all[0,:],transform=crt.crs.PlateCarree(),cmap=cmap, vmin=vrange1[0], vmax=vrange1[1]); 
        bs = ax0.scatter(x=coord_x[st_id],y=coord_y[st_id],s=80,c='red',transform=crt.crs.PlateCarree(),facecolors='none', edgecolors='r'); 
        tmp1 = ax1.scatter([1,2,5,10,25,50,100],ds_gtsm_eva1.loc[(ds_gtsm_eva1['station']==stations[st_id]),['1','2','5','10','25','50','100']],label=f"{settings['yearmin1']}-{settings['yearmax1']}"); ax1.semilogx(); plt.grid();
        tmp2 = ax1.scatter([1,2,5,10,25,50,100],ds_gtsm_eva2.loc[(ds_gtsm_eva1['station']==stations[st_id]),['1','2','5','10','25','50','100']],c='red',label=f"{settings['yearmin2']}-{settings['yearmax2']}");
        
        ax1.legend(loc="upper left")
        ax1.set_xlabel('Return period in years')
        ax1.set_ylabel('Total water level [m]')
        ax1.grid()
        
        ts = ax2.plot(ds_gtsm['time'].sel(stations=stations[st_id]),ds_gtsm['waterlevel']);
        ax2.set_ylabel('Total water level [m]')
        ax0.title.set_text(f"Location {stations[st_id]}")
        ax1.title.set_text(f"Return values")
        figname = f'EVA_station_{str(stations[st_id]).zfill(5)}_pot_plots.png' 
        fig.savefig(f'{dir_eva}/{figname}')
    

#    # Plotting difference of present calculation and dataset on CDS - rv values
#    cmap = mpl.colormaps['magma_r'].resampled(20)
#    vrange=[0,0.5]    
#    print('Plotting... ')
#    fig = plt.figure(figsize=(26,20))
#    axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
#    axs = AxesGrid(fig, 111, axes_class=axes_class,
#            nrows_ncols=(2, 2),
#            share_all=True,
#            axes_pad=1.7,
#            cbar_location='right',
#            cbar_mode='each',
#            cbar_size='3%',
#            cbar_pad=1.0,
#            label_mode='keep')
#    for ii in range(0,4):
#        ax = global_map(axs[ii])
#        bs = ax.scatter(x=coord_x,y=coord_y,
#                        s=25,c=abs(ds_gtsm_eva3[str(rps[ii])].values-ds_rv[ii,:]),transform=crt.crs.PlateCarree(),
#                        cmap=cmap, vmin=vrange[0], vmax=vrange[1])
#        cbar = ax.cax.colorbar(bs)
#    figname = 'EVA_map_compare_dataset_to_CDS_1985-2014_10min.png' 
#    fig.savefig(f'{dir_eva}/{figname}')

    # # Plotting difference of present calculation and dataset on CDS - percentiles
    # cmap = mpl.colormaps['magma_r'].resampled(20)
    # vrange=[0,0.5]    
    # print('Plotting... ')
    # fig = plt.figure(figsize=(26,20))
    # axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    # axs = AxesGrid(fig, 111, axes_class=axes_class,
    #         nrows_ncols=(1, 2),
    #         share_all=True,
    #         axes_pad=1.7,
    #         cbar_location='right',
    #         cbar_mode='each',
    #         cbar_size='3%',
    #         cbar_pad=1.0,
    #         label_mode='keep')
    # ax = global_map(axs[0])
    # bs = ax.scatter(x=coord_x,y=coord_y,
    #                     s=25,c=abs(quan90-ds_perc90['water_level_percentile_90'].values),transform=crt.crs.PlateCarree(),
    #                     cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    # cbar = ax.cax.colorbar(bs)
    # ax = global_map(axs[1])
    # bs = ax.scatter(x=coord_x,y=coord_y,
    #                     s=25,c=abs(quan95-ds_perc95['water_level_percentile_95'].values),transform=crt.crs.PlateCarree(),
    #                     cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    # cbar = ax.cax.colorbar(bs)
    # figname = 'EVA_map_compare_dataset_to_CDS_1985-2014_10min_percentiles.png' 
    # fig.savefig(f'{dir_eva}/{figname}')    

    # # analyze the differences
    # diff = abs(ds_gtsm_eva3['100'].values-ds_rv100['return_mean_water_level'].values)
    # ind = [i for i,v in enumerate(diff) if v > 0.1]

    # data = {'station': ind,
    #         'coord_x': coord_x[ind],
    #         'coord_y': coord_y[ind],
    #         'perc95_new': quan95[ind],
    #         'perc95_old': ds_perc95['water_level_percentile_95'].values[ind],
    #         'rv100 diff': diff[ind],
    #         'rv100 new': ds_gtsm_eva3['100'].values[ind],
    #         'rv100 old': ds_rv100['return_mean_water_level'].values[ind]}
    # df = pd.DataFrame(data)

    # # try to work with timeseries that were not detrended and use numpy for percentiles

    # from eva import detrend, stats
    
    # ds_gtsm0 = xr.open_dataset(f'{dir_eva3}/ds_GTSM-ERA5_1985_2014_stations_00000-00999.nc'); ds.close()

    # # detrend
    # ds_gtsm['sea_level'] = ds_gtsm0['waterlevel']
    # ds_gtsm = ds_gtsm.drop(['waterlevel'])
    # ds_gtsm = detrend(ds_gtsm)
    
    # # compute stats
    # prcts = [0.90,0.95,0.99,0.999]    
    # ds_stats0 = ds_gtsm0.waterlevel.quantile(prcts, dim=('time'))  
    # ds_stats = ds_gtsm.sea_level_detrended.quantile(prcts, dim=('time')) 

    # data = {'station': range(0,1000),
    #         'coord_x': coord_x[0:1000],
    #         'coord_y': coord_y[0:1000],
    #         'perc95_new': quan95[0:1000], 
    #         'perc95_new0': ds_stats0.values[1,:],
    #         'perc95_old': ds_perc95['water_level_percentile_95'].values[0:1000]}
    
    # df2 = pd.DataFrame(data)

    # # Plotting difference of present calculation and dataset on CDS - percentiles - without detrend
    # cmap = mpl.colormaps['hot_r'].resampled(20)
    # vrange=[0,0.1]    
    # print('Plotting... ')
    # fig = plt.figure(figsize=(26,20))
    # axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    # axs = AxesGrid(fig, 111, axes_class=axes_class,
    #         nrows_ncols=(1, 2),
    #         share_all=True,
    #         axes_pad=1.7,
    #         cbar_location='right',
    #         cbar_mode='each',
    #         cbar_size='3%',
    #         cbar_pad=1.0,
    #         label_mode='keep')
    # ax = global_map(axs[0])
    # bs = ax.scatter(x=coord_x[0:1000],y=coord_y[0:1000],
    #                     s=25,c=abs(quan95[0:1000]-ds_perc95['water_level_percentile_95'].values[0:1000]),transform=crt.crs.PlateCarree(),
    #                     cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    # cbar = ax.cax.colorbar(bs)
    # ax = global_map(axs[1])
    # bs = ax.scatter(x=coord_x[0:1000],y=coord_y[0:1000],
    #                     s=25,c=abs(ds_stats0.values[1,:]-ds_perc95['water_level_percentile_95'].values[0:1000]),transform=crt.crs.PlateCarree(),
    #                     cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    # cbar = ax.cax.colorbar(bs)
    # figname = 'EVA_map_compare_dataset_to_CDS_1985-2014_10min_95perc_with_and_without_detrend.png' 
    # fig.savefig(f'{dir_eva}/{figname}')    
