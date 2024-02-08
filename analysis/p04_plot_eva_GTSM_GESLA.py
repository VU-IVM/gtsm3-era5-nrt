# ---
# Author: N. Aleksandrova
# Contact: natalia.aleksandrova@deltares.nl
# Date created: July 2023
# Remarks: gtsm_eva

import pandas as pd
from global_map import global_map
import warnings
warnings.filterwarnings('ignore')
import xarray as xr
import numpy as np
import sys
import os
from math import sqrt, cos, radians
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
sys.path.append("..")
from path_dict import path_dict

from scipy.stats.stats import pearsonr

if __name__ == "__main__":   
    # EVA was made available for three periods:
    settings = {'yearmin1': 1985,'yearmax1': 2014,'mode1':'1hr',
                'yearmin2': 1979,'yearmax2': 2018,'mode2':'1hr',
                'yearmin3': 1950,'yearmax3': 2022,'mode3':'1hr'}
    rps =[1,10,50,100]    
    
    # location of EVA data
    dir_postproc = path_dict['postproc']
    dir_eva_main = os.path.join(dir_postproc,'EVA-GTSM-ERA5')
    dir_eva = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_1950_2022_1hr_v2') 
    
    #locate .csv files
    filenames = 'ds_GTSM-ERA5_1950_2022_stations*eva.csv'
    dir_data = os.path.join(dir_eva,filenames)
    file_list = glob.glob(dir_data)
    file_list.sort()

    #read and merge eva data in one data frame
    ds_gtsm_eva = pd.read_csv(file_list[0])

    for ii in range(1,len(file_list)):
        tmp = pd.read_csv(file_list[ii])
        ds_gtsm_eva = pd.concat([ds_gtsm_eva,tmp])
    del tmp

    #locate and read files that contain coordinates 
    file_list_nc = [str(file_list[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list))]
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
    file_list_nc = [str(file_list[ii]).replace('_eva.csv','_stats.nc') for ii in range(0,len(file_list))]
    ds = xr.open_dataset(file_list_nc[0]); ds.close()
    quan90 = ds['sea_level_detrended'].values[5,:]
    quan95 = ds['sea_level_detrended'].values[6,:]
    for ii in range(1,len(file_list_nc)):
        ds = xr.open_dataset(file_list_nc[ii]); ds.close()
        quan90 = np.append(quan90,ds['sea_level_detrended'].values[5,:])
        quan95 = np.append(quan95,ds['sea_level_detrended'].values[6,:])
    del ds, file_list_nc
    
    # PLOTTING inputs
    cmap = mpl.colormaps['viridis'].resampled(20)
    cmap2 = mpl.colormaps['hot_r']#.resampled(20)
    cmap3 = mpl.colormaps['magma_r'].resampled(20)
    cmap4 = mpl.colormaps['seismic']#.resampled(20)

    # read GESLA dataset
    file_gesla = os.path.join('/gpfs/work1/0/einf3499/','data','ds_gesla_1950_2022_allstations_50yr_max25prt_missing.nc')
    ds_ges = xr.open_dataset(file_gesla); ds_ges.close()       

    # read GTSM dataset coordinates
    
    dir_wlts = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly/waterlevel/'
    dir_wlts2 = '/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/timeseries-GTSM-ERA5-hourly-1979-2018/waterlevel/'

    file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_1950_01_v1.nc')
    ds_gtsm = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds_gtsm.close()
    ds_gtsm.load()

    # find GTSM locations corresponding to GESLA stations
    ids_gtsm=[]
    ids_ges = []
    for ii in range(0,len(ds_ges.stations.values)):
        abslat = np.abs(ds_ges.station_y_coordinate.values[ii] - ds_gtsm.station_y_coordinate.values)
        abslon = np.abs(ds_ges.station_x_coordinate.values[ii] - ds_gtsm.station_x_coordinate.values)
        id_gtsm=np.argmin(abslon**2 + abslat**2)
        x = (radians(ds_ges.station_x_coordinate.values[ii]) - radians(ds_gtsm.station_x_coordinate.values[id_gtsm])) * cos(0.5 * (radians(ds_ges.station_y_coordinate.values[ii]) + radians(ds_gtsm.station_y_coordinate.values[id_gtsm])))
        y = radians(ds_ges.station_y_coordinate.values[ii]) - radians(ds_gtsm.station_y_coordinate.values[id_gtsm])
        d = 6371 * sqrt(x*x + y*y)
        if d < 10: # only consider stations less than 10 km apart
            ids_gtsm.append(id_gtsm)
            ids_ges.append(ii)

    # get a list of all stations in the GTSM dataset
    stations = ds_gtsm.stations.values

    # plot GESLA and GTSM locations
    fig = plt.figure(figsize=(15,10))
    axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(1, 1), label_mode='keep')
    ax = global_map(axs[0])
    #bs = ax.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,color='blue',transform=crt.crs.PlateCarree()); 
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values[ids_gtsm],y=ds_gtsm.station_y_coordinate.values[ids_gtsm],s=15,color='red',transform=crt.crs.PlateCarree()); 
    
    #load GTSM timeseries for selected stations
    # for year in range(1950,2023):
    #     print(f'loading {year}')
    #     for mnth in range(1,13):
    #         print(f'loading {year}, month {mnth}')
    #         if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
    #             file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
    #             ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
    #             if "station_x_coordinate" in list(ds.data_vars):
    #                 ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))                 
    #         else:
    #             file_nc = os.path.join(dir_wlts2,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
    #             ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
    #         ds=ds.sel(stations=stations[ids_gtsm],drop=True)
    #         ds.load()
    #         if ((year == 1950) & (mnth == 1)):
    #             ds_gtsm = ds
    #         else:
    #             ds_gtsm = xr.concat([ds_gtsm,ds],dim="time")
        
    # # save timeseries selection 
    # ds_gtsm.to_netcdf("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1950_2022_1hr_selected_stations_GESLA.nc");

    # open timeseries selection and detrend
    ds_gtsm = xr.open_dataset("/gpfs/work1/0/einf3499/06_model_runs/03_postprocessing/EVA-GTSM-ERA5/period_1950_2022_1hr_selected_stations_GESLA.nc");
    ds_gtsm['sea_level'] = ds_gtsm['waterlevel']
    ds_gtsm = ds_gtsm.drop(['waterlevel'])
    ds_gtsm = detrend(ds_gtsm)
    ds_gtsm = ds_gtsm.drop(['sea_level'])
    ds_gtsm = ds_gtsm.chunk({"time": -1, "stations": "auto"})
    ds_gtsm.load()
    ds_gtsm=ds_gtsm.set_coords(("station_x_coordinate", "station_y_coordinate"))

    # detrend GESLA timeseries
    ds_ges = detrend(ds_ges)
    ds_ges = ds_ges.drop(['sea_level'])
    ds_ges.load()

    vrange1=[-0.5,0.5]; vrange2=[0,0.02]; 

    make_plot = 1

    # create arrays for statistics
    all_rp100 = np.zeros(shape=(len(ids_gtsm), 4))
    all_rp10 = np.zeros(shape=(len(ids_gtsm), 4))
    all_q95 = np.zeros(shape=(len(ids_gtsm), 4))
    all_q99 = np.zeros(shape=(len(ids_gtsm), 4))
    all_ts_rmse = np.zeros(shape=(len(ids_gtsm), 2))
    all_ts_corr = np.zeros(shape=(len(ids_gtsm), 2))
    
    # make overview plots per location
    for ss in range(0,len(ids_gtsm)):
        st_id_gtsm = ids_gtsm[ss]
        st_num_ges = ids_ges[ss]

        print('processing station ',ss,' out of ', len(ids_gtsm))

        ts_gesla = ds_ges.sea_level_detrended.isel(stations=st_num_ges)
        ts_gtsm = ds_gtsm.sea_level_detrended.sel(stations=stations[st_id_gtsm]).sel(time=slice('01-01-1950','31-12-2020'))
        #ts_gesla_ori = ds_ges.sea_level.isel(stations=st_num_ges) - ds_ges.sea_level.isel(stations=st_num_ges).mean() # without detrending, but removing mean
        #ts_gtsm_ori = ds_gtsm.sea_level_detrended.sel(stations=stations[st_id_gtsm]).sel(time=slice('01-01-1950','31-12-2020')) - ds_gtsm.sea_level_detrended.sel(stations=stations[st_id_gtsm]).sel(time=slice('01-01-1950','31-12-2020')).mean()
        
        if np.ndim(ts_gtsm.stations.values)>0:
            ts_gtsm = ts_gtsm.isel(stations=0)
            #ts_gtsm_ori = ts_gtsm_ori.isel(stations=0)
        
        # take only part of the timeseries that overlaps with existing GESLA observations
        mask = 1-(np.isnan(ts_gesla.values) | np.isnan(ts_gtsm.values))
        ts_gtsm = ts_gtsm.where(mask, np.nan)
        ts_gesla = ts_gesla.where(mask, np.nan)
        #ts_gtsm_ori = ts_gtsm_ori.where(mask, np.nan)
        #ts_gesla_ori = ts_gesla_ori.where(mask, np.nan)

        # calculate quantiles
        all_q95[ss] = [float(ts_gtsm.quantile(0.95)), float(ts_gtsm.sel(time=slice('01-01-1979','31-12-2018')).quantile(0.95)), float(ts_gesla.quantile(0.95)), float(ts_gesla.sel(time=slice('01-01-1979','31-12-2018')).quantile(0.95))]
        all_q99[ss] = [float(ts_gtsm.quantile(0.99)), float(ts_gtsm.sel(time=slice('01-01-1979','31-12-2018')).quantile(0.99)), float(ts_gesla.quantile(0.99)), float(ts_gesla.sel(time=slice('01-01-1979','31-12-2018')).quantile(0.99))]

        # calculate timeseries RMSE and correlation
        all_ts_corr[ss] = [pearsonr(ts_gtsm.values[~np.isnan(ts_gtsm.values)], ts_gesla.values[~np.isnan(ts_gesla.values)]).statistic, pearsonr(ts_gtsm.sel(time=slice('01-01-1979','31-12-2018')).values[~np.isnan(ts_gtsm.sel(time=slice('01-01-1979','31-12-2018')).values)], ts_gesla.sel(time=slice('01-01-1979','31-12-2018')).values[~np.isnan(ts_gesla.sel(time=slice('01-01-1979','31-12-2018')).values)]).statistic]
        all_ts_rmse[ss] = [np.sqrt(((ts_gtsm.values[~np.isnan(ts_gtsm.values)] - ts_gesla.values[~np.isnan(ts_gesla.values)]) ** 2).mean()), np.sqrt(((ts_gtsm.sel(time=slice('01-01-1979','31-12-2018')).values[~np.isnan(ts_gtsm.sel(time=slice('01-01-1979','31-12-2018')).values)] - ts_gesla.sel(time=slice('01-01-1979','31-12-2018')).values[~np.isnan(ts_gesla.sel(time=slice('01-01-1979','31-12-2018')).values)]) ** 2).mean())]
        
        # prepare EVA outputs
        # GTSM full timeseries
        var = ts_gtsm.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
        varth = var.quantile(0.99); 
        model_gtsm = EVA(var); del var
        model_gtsm.get_extremes("POT", threshold=varth, r="72H");  del varth
        model_gtsm.fit_model(distribution='genpareto',model='MLE')

        # GESLA full timeseries
        var = ts_gesla.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
        varth = var.quantile(0.99)
        model_gesla = EVA(var); del var
        model_gesla.get_extremes("POT", threshold=varth, r="72H"); del varth
        model_gesla.fit_model(distribution='genpareto',model='MLE')

        # GTSM short timeseries
        var = ts_gtsm.sel(time=slice('01-01-1979','31-12-2018')).to_dataframe().loc[:, 'sea_level_detrended'].dropna()
        varth = var.quantile(0.99); 
        model_gtsm_short = EVA(var); del var
        model_gtsm_short.get_extremes("POT", threshold=varth, r="72H");  del varth
        model_gtsm_short.fit_model(distribution='genpareto',model='MLE')

        # GESLA short timeseries
        var = ts_gesla.sel(time=slice('01-01-1979','31-12-2018')).to_dataframe().loc[:, 'sea_level_detrended'].dropna()
        varth = var.quantile(0.99)
        model_gesla_short = EVA(var); del var
        model_gesla_short.get_extremes("POT", threshold=varth, r="72H"); del varth
        model_gesla_short.fit_model(distribution='genpareto',model='MLE')

        # get pandas with EVA results for plotting
        model_summary_gtsm = model_gtsm.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)
        model_summary_gesla = model_gesla.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)
        model_summary_gtsm_short = model_gtsm_short.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)
        model_summary_gesla_short = model_gesla_short.get_summary(return_period=[1,2,5,10,50,100,500,1000],alpha=0.95)

        all_rp10[ss] = [model_summary_gtsm['return value'].loc[10],model_summary_gtsm_short['return value'].loc[10], model_summary_gesla['return value'].loc[10], model_summary_gesla_short['return value'].loc[10]]

        all_rp100[ss] = [model_summary_gtsm['return value'].loc[100],model_summary_gtsm_short['return value'].loc[100], model_summary_gesla['return value'].loc[100], model_summary_gesla_short['return value'].loc[100]]
        
        if make_plot:
            # plot
            fig = plt.figure(figsize=(20,20))
            ax0 = plt.subplot2grid((4, 2), (0, 0), colspan=1, projection=crt.crs.Robinson())
            ax1 = plt.subplot2grid((4, 2), (0, 1), colspan=1, rowspan=1)
            ax2 = plt.subplot2grid((4, 2), (1, 0), colspan=1, rowspan=1)
            ax3 = plt.subplot2grid((4, 2), (1, 1), colspan=1, rowspan=1,sharex=ax2,sharey=ax2)
            ax4 = plt.subplot2grid((4, 2), (2, 0), colspan=2, rowspan=1)
            ax5 = plt.subplot2grid((4, 2), (3, 0), colspan=1, rowspan=1)
            ax6 = plt.subplot2grid((4, 2), (3, 1), colspan=1, rowspan=1, projection=crt.crs.PlateCarree())

            #    plot location on the map
            ax0 = global_map(ax0)
            bs0 = ax0.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,color='blue',transform=crt.crs.PlateCarree()); # all considered GESLA stations
            bs = ax0.scatter(x=ts_gtsm.station_x_coordinate.values,y=ts_gtsm.station_y_coordinate.values,marker ='.',s=200,transform=crt.crs.PlateCarree(),facecolors='none',edgecolors='red',linewidth=3); 
            ax0.title.set_text(f"Location {int(stations[st_id_gtsm])}")
        
            #   plot return values on a log plot per dataset
            model_gtsm.plot_return_values(return_period=[1, 2, 5, 10, 25, 50, 100,500,1000],alpha=0.95,ax=ax2)
            model_gesla.plot_return_values(return_period=[1, 2, 5, 10, 25, 50, 100,500,1000],alpha=0.95,ax=ax3)
            ax2.set_ylabel('Still water level [m]'); ax3.set_ylabel('Still water level [m]')
            ax2.title.set_text(f"Return values based on GTSM-ERA5 dataset for 1950-2022"); ax3.title.set_text(f"Return values based on GESLA dataset"); 

            # plot EVA confidence intervals together for comparison
            ax1.semilogx()
            ax1.grid(True, which="both")
            ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))
        
            for col in ["lower ci", "upper ci"]:
                ax1.plot(model_summary_gtsm.index.values,model_summary_gtsm.loc[:, col].values,color="#5199FF",lw=1,ls="--",zorder=14)
                ax1.plot(model_summary_gesla.index.values,model_summary_gesla.loc[:, col].values,color="#ff6e51",lw=1,ls="--",zorder=15)
            ax1.fill_between(model_summary_gtsm.index.values,model_summary_gtsm.loc[:, "lower ci"].values,model_summary_gtsm.loc[:, "upper ci"].values, facecolor="#5199FF",edgecolor="None",alpha=0.25,zorder=9)
            ax1.fill_between(model_summary_gesla.index.values,model_summary_gesla.loc[:, "lower ci"].values,model_summary_gesla.loc[:, "upper ci"].values, facecolor="#ff6e51",edgecolor="None",alpha=0.25,zorder=10)
            ax1.plot(model_summary_gtsm.index.values,model_summary_gtsm.loc[:, "return value"].values,color="#285fad",lw=2,ls="-",zorder=24,label='Based on GTSM-ERA5 dataset')
            ax1.plot(model_summary_gesla.index.values,model_summary_gesla.loc[:, "return value"].values,color="#b8391f",lw=2,ls="-",zorder=25,label='Based on GESLA dataset')
            ax1.set_xlabel("Return period")
            ax1.set_ylabel("Still water level [m]")
            ax1.legend(loc="upper left")
            ax1.title.set_text(f"Effect of timeseries length on modelled return values")

            # plot timeseries
            ts = ax4.plot(ts_gesla['time'].values,ts_gesla.values,'r.',label='GESLA');
            ts = ax4.plot(ts_gtsm['time'].values,ts_gtsm.values,'b-',alpha=0.5,label='GTSM-ERA5');
            ax4.set_ylabel('Still water level [m]')
            ax4.grid(); ax4.title.set_text('Full timeseries')

            # plot density scatter
            ax5.hist2d(ts_gtsm.values[~np.isnan(ts_gtsm.values)], ts_gesla.values[~np.isnan(ts_gesla.values)], bins=(300, 30), cmap='Blues')
            tmp = np.nanmax(ts_gesla.values)
            ax5.set(xlim=(-tmp, tmp), ylim=(-tmp, tmp))
            ax5.plot([-tmp,tmp],[-tmp,tmp])
            ax5.set_aspect('equal')  
            corr = pearsonr(ts_gtsm.values[~np.isnan(ts_gtsm.values)], ts_gesla.values[~np.isnan(ts_gesla.values)])
            ax5.text(tmp,-tmp,f'corr={np.round(corr.statistic,3)}', verticalalignment ='bottom', horizontalalignment ='right', fontsize = 14, fontweight ='bold'); 

            bs0 = ax6.scatter(x=ts_gesla.station_x_coordinate.values,y=ts_gesla.station_y_coordinate.values,s=15,color='blue',transform=crt.crs.PlateCarree()); 
            bs = ax6.scatter(x=ts_gtsm.station_x_coordinate.values,y=ts_gtsm.station_y_coordinate.values,s=15,color='red',transform=crt.crs.PlateCarree()); 
            ax6.coastlines()
            ax6.set_extent([ts_gtsm.station_x_coordinate.values-0.1, ts_gtsm.station_x_coordinate.values+0.1, ts_gtsm.station_y_coordinate.values-0.1, ts_gtsm.station_y_coordinate.values+0.1])
            ax6.title.set_text(f"Location map")

            del ts_gtsm, ts_gesla, model_summary_gtsm, model_summary_gesla

            figname = f'EVA_station_{str(int(stations[st_id_gtsm])).zfill(5)}_RV_1950-2022_vs_GESLA.png' 
            fig.savefig(f'{dir_eva_main}/timeseries_plots_GESLA/{figname}')
            plt.clf()
    

    # print statistics

    print(f'Hourly timeseries: RMSE, GTSM-ERA5-E: {np.mean(all_ts_rmse[:,0]):.3} ({np.std(all_ts_rmse[:,0]):.3})')
    print(f'Hourly timeseries: RMSE, GTSM-ERA5: {np.mean(all_ts_rmse[:,1]):.3} ({np.std(all_ts_rmse[:,1]):.3})')

    print(f'Hourly timeseries: pearson nr, GTSM-ERA5-E: {np.mean(all_ts_corr[:,0]):.3} ({np.std(all_ts_corr[:,0]):.3})')
    print(f'Hourly timeseries: pearson nr, GTSM-ERA5: {np.mean(all_ts_corr[:,1]):.3} ({np.std(all_ts_corr[:,1]):.3})')
    
    print('  ')
    print(f'95th percentile: Mean bias GTSM-ERA5-E: {np.mean(all_q95[:,0]-all_q95[:,2]):.3} ({np.std(all_q95[:,0]-all_q95[:,2]):.3})')
    print(f'95th percentile: Mean bias GTSM-ERA5: {np.mean(all_q95[:,1]-all_q95[:,3]):.3} ({np.std(all_q95[:,1]-all_q95[:,3]):.3})')
    print(f'95th percentile: Mean abs error GTSM-ERA5-E: {np.mean(abs(all_q95[:,0]-all_q95[:,2])):.3} ({np.std(abs(all_q95[:,0]-all_q95[:,2])):.3})')
    print(f'95th percentile: Mean abs error GTSM-ERA5: {np.mean(abs(all_q95[:,1]-all_q95[:,3])):.3} ({np.std(abs(all_q95[:,1]-all_q95[:,3])):.3})')
    print(f'95th percentile: Mean abs perc error GTSM-ERA5-E: {np.mean(abs(all_q95[:,0]-all_q95[:,2])/all_q95[:,2])*100:.3} ({np.std(abs(all_q95[:,0]-all_q95[:,2])/all_q95[:,2])*100:.3})')
    print(f'95th percentile: Mean abs perc error GTSM-ERA5: {np.mean(abs(all_q95[:,1]-all_q95[:,3])/all_q95[:,3]):.3} ({np.std(abs(all_q95[:,1]-all_q95[:,3])/all_q95[:,3]):.3})')

    print('  ')
    print(f'99th percentile: Mean bias GTSM-ERA5-E: {np.mean(all_q99[:,0]-all_q99[:,2]):.3} ({np.std(all_q99[:,0]-all_q99[:,2]):.3})')
    print(f'99th percentile: Mean bias GTSM-ERA5: {np.mean(all_q99[:,1]-all_q99[:,3]):.3} ({np.std(all_q99[:,1]-all_q99[:,3]):.3})')
    print(f'99th percentile: Mean abs error GTSM-ERA5-E: {np.mean(abs(all_q99[:,0]-all_q99[:,2])):.3} ({np.std(abs(all_q99[:,0]-all_q99[:,2])):.3})')
    print(f'99th percentile: Mean abs error GTSM-ERA5: {np.mean(abs(all_q99[:,1]-all_q99[:,3])):.3} ({np.std(abs(all_q99[:,1]-all_q99[:,3])):.3})')
    print(f'99th percentile: Mean abs perc error GTSM-ERA5-E: {np.mean(abs(all_q99[:,0]-all_q99[:,2])/all_q99[:,2])*100:.3} ({np.std(abs(all_q99[:,0]-all_q99[:,2])/all_q99[:,2])*100:.3})')
    print(f'99th percentile: Mean abs perc error GTSM-ERA5: {np.mean(abs(all_q99[:,1]-all_q99[:,3])/all_q99[:,3]):.3} ({np.std(abs(all_q99[:,1]-all_q99[:,3])/all_q99[:,3]):.3})')

    print('  ')
    print(f'RP10: Mean bias GTSM-ERA5-E: {np.mean(all_rp10[:,0]-all_rp10[:,2]):.3} ({np.std(all_rp10[:,0]-all_rp10[:,2]):.3})')
    print(f'RP10: Mean bias GTSM-ERA5: {np.mean(all_rp10[:,1]-all_rp10[:,3]):.3} ({np.std(all_rp10[:,1]-all_rp10[:,3]):.3})')
    print(f'RP10: Mean abs error GTSM-ERA5-E: {np.mean(abs(all_rp10[:,0]-all_rp10[:,2])):.3} ({np.std(abs(all_rp10[:,0]-all_rp10[:,2])):.3})')
    print(f'RP10: Mean abs error GTSM-ERA5: {np.mean(abs(all_rp10[:,1]-all_rp10[:,3])):.3} ({np.std(abs(all_rp10[:,1]-all_rp10[:,3])):.3})')
    print(f'RP10: Mean abs perc error GTSM-ERA5-E: {np.mean(abs(all_rp10[:,0]-all_rp10[:,2])/all_rp10[:,2])*100:.3} ({np.std(abs(all_rp10[:,0]-all_rp10[:,2])/all_rp10[:,2])*100:.3})')
    print(f'RP10: Mean abs perc error GTSM-ERA5: {np.mean(abs(all_rp10[:,1]-all_rp10[:,3])/all_rp10[:,3]):.3} ({np.std(abs(all_rp10[:,1]-all_rp10[:,3])/all_rp10[:,3]):.3})')

    print('  ')
    print(f'RP100: Mean bias GTSM-ERA5-E: {np.mean(all_rp100[:,0]-all_rp100[:,2]):.3} ({np.std(all_rp100[:,0]-all_rp100[:,2]):.3})')
    print(f'RP100: Mean bias GTSM-ERA5: {np.mean(all_rp100[:,1]-all_rp100[:,3]):.3} ({np.std(all_rp100[:,1]-all_rp100[:,3]):.3})')
    print(f'RP100: Mean abs error GTSM-ERA5-E: {np.mean(abs(all_rp100[:,0]-all_rp100[:,2])):.3} ({np.std(abs(all_rp100[:,0]-all_rp100[:,2])):.3})')
    print(f'RP100: Mean abs error GTSM-ERA5: {np.mean(abs(all_rp100[:,1]-all_rp100[:,3])):.3} ({np.std(abs(all_rp100[:,1]-all_rp100[:,3])):.3})')
    print(f'RP100: Mean abs perc error GTSM-ERA5-E: {np.mean(abs(all_rp100[:,0]-all_rp100[:,2])/all_rp100[:,2])*100:.3} ({np.std(abs(all_rp100[:,0]-all_rp100[:,2])/all_rp100[:,2])*100:.3})')
    print(f'RP100: Mean abs perc error GTSM-ERA5: {np.mean(abs(all_rp100[:,1]-all_rp100[:,3])/all_rp100[:,3]):.3} ({np.std(abs(all_rp100[:,1]-all_rp100[:,3])/all_rp100[:,3]):.3})')


    
    # make plot overview for q95 and RP100

    fig = plt.figure(figsize=(12,10))
    axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(2, 2), share_all=True, axes_pad=1,cbar_location='right',cbar_mode='each',cbar_size='3%',cbar_pad=0.3, label_mode='keep')
    
    # plot q95 and RP100 map
    ax = global_map(axs[0])
    bs = ax.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,c=all_q95[:,0],cmap=cmap,transform=crt.crs.PlateCarree()); #vmin=vrange1[0], vmax=vrange1[1]
    cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=14)
    ax.set_title('95th percentile values based on GTSM-ERA5-E')

    ax = global_map(axs[1])
    bs = ax.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,c=all_rp100[:,0],cmap=cmap3,transform=crt.crs.PlateCarree()); #vmin=vrange1[0], vmax=vrange1[1]
    cbar = ax.cax.colorbar(bs); cbar.set_label('Still water level [m]',fontsize=14)
    ax.set_title('100-year return values based on GTSM-ERA5-E')
    
    # plot q95 and RP100 bias
    ax = global_map(axs[2])
    bs =ax.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,c=all_q95[:,0]-all_q95[:,2],cmap=cmap4,transform=crt.crs.PlateCarree(),vmin=-1, vmax=1);
    cbar = ax.cax.colorbar(bs); cbar.set_label('Bias [m]',fontsize=14)
    ax.set_title('Difference in 95th percentile values \n for GTSM-ERA5-E vs. observations')

    ax = global_map(axs[3])
    bs =ax.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,c=all_rp100[:,0]-all_rp100[:,2],cmap=cmap4,transform=crt.crs.PlateCarree(),vmin=-1, vmax=1); #
    cbar = ax.cax.colorbar(bs); cbar.set_label('Bias [m]',fontsize=14)
    ax.set_title('Difference in 100-year return values \n for GTSM-ERA5-E vs. observations')


