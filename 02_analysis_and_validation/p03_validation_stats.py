# ---
# Author: N. Aleksandrova
# Contact: natalia.aleksandrova@deltares.nl
# Date created: July 2023
# Remarks: gtsm_eva
#%%
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
from pyextremes import EVA
#from path_dict import path_dict
from scipy.stats.stats import pearsonr
sys.path.append("..")

#%%
if __name__ == "__main__":   
 
    dir_eva_main = r'p:\archivedprojects\11210221-gtsm-reanalysis\GTSM-ERA5-E_dataset\EVA-GTSM-ERA5'

    dir_out = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\validation_plots'
    dir_out_stats = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\validation_stats'
    
    # Get GESLA data
    filename_geslaselection = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\period_1950_2024_1hr_selected_stations_GESLA_filtered.nc'
    stationlist = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\overview\waterlevel_data_netcdf_overview_filtered_75perc_bothperiods.csv'
    geslalist=pd.read_csv(stationlist)

    # read GESLA dataset
    #file_gesla = os.path.join('/gpfs/work1/0/einf3499/','data','ds_gesla_1950_2022_allstations_50yr_max25prt_missing.nc')
    #file_gesla = r'c:\Users\aleksand\OneDrive - Stichting Deltares\Documents\Projects\GTSM-ERA5\paper_GTSM-ERA5-E\analysis\gesla\ds_gesla_1950_2022_allstations_50yr_max25prt_missing.nc'
    file_gesla = r'p:\archivedprojects\11210221-gtsm-reanalysis\GESLA\ds_gesla_1950_2022_allstations_50yr_max25prt_missing.nc'
    ds_ges = xr.open_dataset(file_gesla);    
    
    # Keep only preselected stations with long records and good coverage
    ds_ges = ds_ges.where(ds_ges.station_name.isin(list(geslalist.station_name)), drop=True)

    # keep only unique stations wih longest records
    non_nan_count = ds_ges.sea_level.notnull().sum(dim="time") 
    df = pd.DataFrame({
        "station": ds_ges.stations.values,
        "station_name": ds_ges.station_name.values,
        "n_valid": non_nan_count.values,
    })
    best_stations = (
        df.sort_values("n_valid", ascending=False)
        .drop_duplicates(subset="station_name")
        ["station"]
    )
    ds_ges = ds_ges.sel(stations=best_stations.values)
    ds_ges = ds_ges.rename({"sea_level": "waterlevel"})
    # read GTSM dataset coordinates
    dir_wlts = r'p:\archivedprojects\11210221-gtsm-reanalysis\GTSM-ERA5-E_dataset\waterlevel'

    file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_1950_01_v3.nc')
    ds_gtsm = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds_gtsm.close()
    ds_gtsm.station_x_coordinate.load()
    ds_gtsm.station_y_coordinate.load()

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

    if not os.path.isfile(filename_geslaselection):
        # plot GESLA and GTSM locations
        fig = plt.figure(figsize=(15,10))
        axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
        axs = AxesGrid(fig, 111, axes_class=axes_class, nrows_ncols=(1, 1), label_mode='keep')
        ax = global_map(axs[0])
        #bs = ax.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,color='blue',transform=crt.crs.PlateCarree()); 
        bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values[ids_gtsm],y=ds_gtsm.station_y_coordinate.values[ids_gtsm],s=15,color='red',transform=crt.crs.PlateCarree()); 
        fig.savefig(os.path.join(dir_eva_main, 'GTSM_stations_close_to_GESLA_stations_alltime.png'))
        
        # check coverage
        numel_short = np.count_nonzero(~np.isnan(ds_ges.isel(stations=ids_ges).sel(time=slice('1979-01-01','2018-01-01'))['sea_level'].values))
        numel_long = np.count_nonzero(~np.isnan(ds_ges.isel(stations=ids_ges).sel(time=slice('1950-01-01','2025-01-01'))['sea_level'].values))
        numel_extra = np.count_nonzero(~np.isnan(ds_ges.isel(stations=ids_ges).sel(time=slice('1950-01-01','1979-01-01'))['sea_level'].values))

        data_coverage_short = numel_short / np.size(ds_ges.isel(stations=ids_ges).sel(time=slice('1979-01-01','2018-01-01'))['sea_level'])*100
        data_coverage_long = numel_long / np.size(ds_ges.isel(stations=ids_ges).sel(time=slice('1950-01-01','2025-01-01'))['sea_level'])*100
        data_coverage_extra = numel_extra / np.size(ds_ges.isel(stations=ids_ges).sel(time=slice('1950-01-01','1979-01-01'))['sea_level'])*100
        
        print(f'Data coverage in the period 1979-2018: {data_coverage_short}')
        print(f'Data coverage in the period 1950-2024: {data_coverage_long}')
        print(f'Data coverage in the period 1950-1978: {data_coverage_extra}')
    
        #load GTSM timeseries for selected stations
        for year in range(1950,2025):
            print(f'loading {year}')
            for mnth in range(1,13):
                print(f'loading {year}, month {mnth}')
                file_nc = os.path.join(dir_wlts,f'reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v3.nc')
                ds = xr.open_dataset(file_nc,chunks={'stations': 1000});
                if "station_x_coordinate" in list(ds.data_vars):
                    ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))                 
                ds=ds.sel(stations=stations[ids_gtsm],drop=True)
                #ds.load()
                if ((year == 1950) & (mnth == 1)):
                    ds_gtsm = ds
                else:
                    ds_gtsm = xr.concat([ds_gtsm,ds],dim="time")
                ds.close(); del ds
        
        # # save timeseries selection 
        filename_geslaselection = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\period_1950_2024_1hr_selected_stations_GESLA_filtered_v2.nc'
        ds_gtsm.to_netcdf(filename_geslaselection)
    else:
        ds_gtsm = xr.open_dataset(filename_geslaselection,chunks={'stations': 30})
   
   #%%
    # Calculate metrics

    make_plot = 1

    # create arrays for statistics
    all_bias_day = np.zeros(shape=(len(ids_gtsm), 2))
    all_mae_day = np.zeros(shape=(len(ids_gtsm), 2))
    all_mape_day = np.zeros(shape=(len(ids_gtsm), 2))
    all_rmse_day = np.zeros(shape=(len(ids_gtsm), 2))
    all_corr_day = np.zeros(shape=(len(ids_gtsm), 2))

    all_bias_mon = np.zeros(shape=(len(ids_gtsm), 2))
    all_mae_mon = np.zeros(shape=(len(ids_gtsm), 2))
    all_mape_mon = np.zeros(shape=(len(ids_gtsm), 2))
    all_rmse_mon = np.zeros(shape=(len(ids_gtsm), 2))
    all_corr_mon = np.zeros(shape=(len(ids_gtsm), 2))    

    all_bias_yr = np.zeros(shape=(len(ids_gtsm), 2))
    all_mae_yr = np.zeros(shape=(len(ids_gtsm), 2))
    all_mape_yr = np.zeros(shape=(len(ids_gtsm), 2))
    all_rmse_yr = np.zeros(shape=(len(ids_gtsm), 2))
    all_corr_yr = np.zeros(shape=(len(ids_gtsm), 2))

    def detrend(ds: xr.DataArray, plot = False):
        ''' remove annual means'''
        ds = ds.assign_coords(year=ds.time.dt.strftime("%Y"))
        ds_new = (ds.groupby("year") - ds.groupby("year").mean("time"))
        ds['waterlevel'] = ds_new['waterlevel'] - ds_new['waterlevel'].mean()
        return ds
    
    for ss in range(0,len(ids_gtsm)):
        st_gtsm = stations[ids_gtsm[ss]]
        st_num_ges = ids_ges[ss]

        print('processing station ',ss,' out of ', len(ids_gtsm))

        ds_gesla_sel = ds_ges.isel(stations=st_num_ges)
        # reference level correction
        ds_gesla_sel['waterlevel'] = ds_gesla_sel['waterlevel'] - np.nanmean(ds_gesla_sel['waterlevel'].sel(time=slice('1986-01-01','2005-12-31')).values) 

        ds_gtsm_sel = ds_gtsm.sel(stations=st_gtsm).sel(time=slice('01-01-1950','31-12-2020'))
        ds_gtsm_sel.load()

        ds_gesla_sel = detrend(ds_gesla_sel)
        ds_gtsm_sel = detrend(ds_gtsm_sel)

        ts_gesla = ds_gesla_sel.waterlevel
        ts_gtsm = ds_gtsm_sel.waterlevel
        del ds_gesla_sel, ds_gtsm_sel

        if np.ndim(ts_gtsm.stations.values)>0:
            ts_gtsm = ts_gtsm.isel(stations=0)

        # take only part of the timeseries that overlaps with existing GESLA observations
        mask = 1-(np.isnan(ts_gtsm.values) | np.isnan(ts_gesla.values))
        ts_gtsm = ts_gtsm.where(mask, np.nan)
        ts_gesla = ts_gesla.where(mask, np.nan)

        # Daily maxima
        ts_gtsm_day = ts_gtsm.resample(time='1D').max()
        ts_gesla_day = ts_gesla.resample(time='1D').max()
        all_corr_day[ss] = [pearsonr(ts_gtsm_day.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gtsm_day.sel(time=slice('01-01-1950','31-12-1978')).values)], ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values)]).statistic, 
                               pearsonr(ts_gtsm_day.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gtsm_day.sel(time=slice('01-01-1979','31-12-2020')).values)], ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values)]).statistic]
        all_rmse_day[ss] = [np.sqrt(((ts_gtsm_day.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gtsm_day.sel(time=slice('01-01-1950','31-12-1978')).values)] - ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values)]) ** 2).mean()), 
                               np.sqrt(((ts_gtsm_day.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gtsm_day.sel(time=slice('01-01-1979','31-12-2020')).values)] - ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values)]) ** 2).mean())]
        all_bias_day[ss] = [np.nanmean(ts_gtsm_day.sel(time=slice('01-01-1950','31-12-1978')).values - ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values),
                              np.nanmean(ts_gtsm_day.sel(time=slice('01-01-1979','31-12-2020')).values - ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values)]
        all_mae_day[ss] = [np.nanmean(np.abs(ts_gtsm_day.sel(time=slice('01-01-1950','31-12-1978')).values - ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values)),
                             np.nanmean(np.abs(ts_gtsm_day.sel(time=slice('01-01-1979','31-12-2020')).values - ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values))]
        all_mape_day[ss] = [np.nanmean(np.abs((ts_gtsm_day.sel(time=slice('01-01-1950','31-12-1978')).values - ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values) / ts_gesla_day.sel(time=slice('01-01-1950','31-12-1978')).values)),
                             np.nanmean(np.abs((ts_gtsm_day.sel(time=slice('01-01-1979','31-12-2020')).values - ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values) / ts_gesla_day.sel(time=slice('01-01-1979','31-12-2020')).values))]

        del ts_gtsm_day, ts_gesla_day

        # Monthly maxima
        ts_gtsm_month = ts_gtsm.resample(time='1ME').max()
        ts_gesla_month = ts_gesla.resample(time='1ME').max()
        all_corr_mon[ss] = [pearsonr(ts_gtsm_month.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gtsm_month.sel(time=slice('01-01-1950','31-12-1978')).values)], ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values)]).statistic, 
                                 pearsonr(ts_gtsm_month.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gtsm_month.sel(time=slice('01-01-1979','31-12-2020')).values)], ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values)]).statistic]
        all_rmse_mon[ss] = [np.sqrt(((ts_gtsm_month.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gtsm_month.sel(time=slice('01-01-1950','31-12-1978')).values)] - ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values[~np.isnan(ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values)]) ** 2).mean()), 
                                 np.sqrt(((ts_gtsm_month.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gtsm_month.sel(time=slice('01-01-1979','31-12-2020')).values)] - ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values[~np.isnan(ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values)]) ** 2).mean())]        
        all_bias_mon[ss] = [np.nanmean(ts_gtsm_month.sel(time=slice('01-01-1950','31-12-1978')).values - ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values), 
                                np.nanmean(ts_gtsm_month.sel(time=slice('01-01-1979','31-12-2020')).values - ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values)]
        all_mae_mon[ss] = [np.nanmean(np.abs(ts_gtsm_month.sel(time=slice('01-01-1950','31-12-1978')).values - ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values)),      
                                np.nanmean(np.abs(ts_gtsm_month.sel(time=slice('01-01-1979','31-12-2020')).values - ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values))]
        all_mape_mon[ss] = [np.nanmean(np.abs((ts_gtsm_month.sel(time=slice('01-01-1950','31-12-1978')).values - ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values) / ts_gesla_month.sel(time=slice('01-01-1950','31-12-1978')).values)),
                                np.nanmean(np.abs((ts_gtsm_month.sel(time=slice('01-01-1979','31-12-2020')).values - ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values) / ts_gesla_month.sel(time=slice('01-01-1979','31-12-2020')).values))]

        # Annual maxima
        ts_gtsm_am = ts_gtsm.groupby(ts_gtsm['time'].dt.year).max()
        ts_gesla_am = ts_gesla.groupby(ts_gesla['time'].dt.year).max()
        all_corr_yr[ss] = [pearsonr(ts_gtsm_am.sel(year=slice(1950,1978)).values[~np.isnan(ts_gtsm_am.sel(year=slice(1950,1978)).values)], ts_gesla_am.sel(year=slice(1950,1978)).values[~np.isnan(ts_gesla_am.sel(year=slice(1950,1978)).values)]).statistic, 
                              pearsonr(ts_gtsm_am.sel(year=slice(1979,2020)).values[~np.isnan(ts_gtsm_am.sel(year=slice(1979,2020)).values)], ts_gesla_am.sel(year=slice(1979,2020)).values[~np.isnan(ts_gesla_am.sel(year=slice(1979,2020)).values)]).statistic]
        all_rmse_yr[ss] = [np.sqrt(((ts_gtsm_am.sel(year=slice(1950,1978)).values[~np.isnan(ts_gtsm_am.sel(year=slice(1950,1978)).values)] - ts_gesla_am.sel(year=slice(1950,1978)).values[~np.isnan(ts_gesla_am.sel(year=slice(1950,1978)).values)]) ** 2).mean()), 
                              np.sqrt(((ts_gtsm_am.sel(year=slice(1979,2020)).values[~np.isnan(ts_gtsm_am.sel(year=slice(1979,2020)).values)] - ts_gesla_am.sel(year=slice(1979,2020)).values[~np.isnan(ts_gesla_am.sel(year=slice(1979,2020)).values)]) ** 2).mean())]

        all_bias_yr[ss] = [np.nanmean(ts_gtsm_am.sel(year=slice(1950,1978)).values - ts_gesla_am.sel(year=slice(1950,1978)).values),
                              np.nanmean(ts_gtsm_am.sel(year=slice(1979,2020)).values - ts_gesla_am.sel(year=slice(1979,2020)).values)] 
        all_mae_yr[ss] = [np.nanmean(np.abs(ts_gtsm_am.sel(year=slice(1950,1978)).values - ts_gesla_am.sel(year=slice(1950,1978)).values)),
                             np.nanmean(np.abs(ts_gtsm_am.sel(year=slice(1979,2020)).values - ts_gesla_am.sel(year=slice(1979,2020)).values))]
        all_mape_yr[ss] = [np.nanmean(np.abs((ts_gtsm_am.sel(year=slice(1950,1978)).values - ts_gesla_am.sel(year=slice(1950,1978)).values) / ts_gesla_am.sel(year=slice(1950,1978)).values)),
                             np.nanmean(np.abs((ts_gtsm_am.sel(year=slice(1979,2020)).values - ts_gesla_am.sel(year=slice(1979,2020)).values) / ts_gesla_am.sel(year=slice(1979,2020)).values))]

        # # Plot AM validation
        # fig,ax=plt.subplots(figsize=(15,7))
        # ts_gesla_am.plot(ax=ax, label='GESLA')
        # ts_gtsm_am.plot(ax=ax, label='GTSM-ERA5')
        # ax.legend()
        # ax.grid()
        # figname = f'station_{str(int(stations[st_id_gtsm])).zfill(5)}_AnnualMaxima_GTSM_vs_GESLA.png' 
        # fig.savefig(f'{dir_eva_main}/timeseries_plots_GESLA/AnnualMaxima/{figname}')
        # plt.clf()
        # del ts_gtsm_am, ts_gesla_am

        if make_plot:
            # plot
            fig = plt.figure(figsize=(20,16))
            ax0 = plt.subplot2grid((3, 2), (0, 0), colspan=1, projection=crt.crs.Robinson())
            ax1 = plt.subplot2grid((3, 2), (0, 1), colspan=1, rowspan=1, projection=crt.crs.PlateCarree())
            ax2 = plt.subplot2grid((3, 2), (1, 0), colspan=2, rowspan=1)
            ax3 = plt.subplot2grid((3, 2), (2, 0), colspan=1, rowspan=1)
            ax4 = plt.subplot2grid((3, 2), (2, 1), colspan=1, rowspan=1)

            #    plot location on the map
            ax0 = global_map(ax0)
            bs0 = ax0.scatter(x=ds_ges.station_x_coordinate.values[ids_ges],y=ds_ges.station_y_coordinate.values[ids_ges],s=15,color='blue',transform=crt.crs.PlateCarree()); # all considered GESLA stations
            bs = ax0.scatter(x=ts_gtsm.station_x_coordinate.values,y=ts_gtsm.station_y_coordinate.values,marker ='.',s=200,transform=crt.crs.PlateCarree(),facecolors='none',edgecolors='red',linewidth=3); 
            ax0.title.set_text(f"Location {int(ts_gtsm.stations.values)}: {ds_ges.station_name.values[ids_ges[ss]]}")
        
            # plot zoomed location
            bs0 = ax1.scatter(x=ts_gesla.station_x_coordinate.values,y=ts_gesla.station_y_coordinate.values,s=15,color='blue',transform=crt.crs.PlateCarree()); 
            bs = ax1.scatter(x=ts_gtsm.station_x_coordinate.values,y=ts_gtsm.station_y_coordinate.values,s=15,color='red',transform=crt.crs.PlateCarree()); 
            ax1.coastlines()
            ax1.stock_img()
            ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            ax1.set_extent([ts_gtsm.station_x_coordinate.values-0.2, ts_gtsm.station_x_coordinate.values+0.2, ts_gtsm.station_y_coordinate.values-0.2, ts_gtsm.station_y_coordinate.values+0.2])
            ax1.title.set_text(f"Location map")

            # plot timeseries
            ts = ax2.plot(ts_gesla['time'].values,ts_gesla.values,color='grey',label='GESLA');
            ts = ax2.plot(ts_gtsm['time'].values,ts_gtsm.values,'b-',alpha=0.5,label='GTSM-ERA5');
            ax2.set_ylabel('Still water level [m]')
            ax2.legend()
            ax2.grid(); ax2.title.set_text('Full timeseries')

            # plot density scatter
            ax3.hist2d(ts_gtsm.sel(time=slice("1950-01-01","1978-12-31")).values[~np.isnan(ts_gtsm.sel(time=slice("1950-01-01","1978-12-31")).values)], 
                       ts_gesla.sel(time=slice("1950-01-01","1978-12-31")).values[~np.isnan(ts_gesla.sel(time=slice("1950-01-01","1978-12-31")).values)], bins=(300, 30), cmap='Blues')
            tmp = np.nanmax(ts_gesla.sel(time=slice("1950-01-01","1978-12-31")).values)
            ax3.set(xlim=(-tmp, tmp), ylim=(-tmp, tmp))
            ax3.plot([-tmp,tmp],[-tmp,tmp])
            ax3.set_aspect('equal')  
            corr = pearsonr(ts_gtsm.sel(time=slice("1950-01-01","1978-12-31")).values[~np.isnan(ts_gtsm.sel(time=slice("1950-01-01","1978-12-31")).values)], ts_gesla.sel(time=slice("1950-01-01","1978-12-31")).values[~np.isnan(ts_gesla.sel(time=slice("1950-01-01","1978-12-31")).values)])
            ax3.text(tmp,-tmp,f'corr={np.round(corr.statistic,3)}', verticalalignment ='bottom', horizontalalignment ='right', fontsize = 14, fontweight ='bold'); 
            ax3.set_xlabel('GTSM-ERA5 [m]'); ax3.set_ylabel('GESLA [m]'); ax3.set_title('1950-1978')

            ax4.hist2d(ts_gtsm.sel(time=slice("1979-01-01","2020-12-31")).values[~np.isnan(ts_gtsm.sel(time=slice("1979-01-01","2020-12-31")).values)], 
                       ts_gesla.sel(time=slice("1979-01-01","2020-12-31")).values[~np.isnan(ts_gesla.sel(time=slice("1979-01-01","2020-12-31")).values)], bins=(300, 30), cmap='Blues')
            tmp = np.nanmax(ts_gesla.sel(time=slice("1979-01-01","2020-12-31")).values)
            ax4.set(xlim=(-tmp, tmp), ylim=(-tmp, tmp)) 
            ax4.plot([-tmp,tmp],[-tmp,tmp])
            ax4.set_aspect('equal')
            corr = pearsonr(ts_gtsm.sel(time=slice("1979-01-01","2020-12-31")).values[~np.isnan(ts_gtsm.sel(time=slice("1979-01-01","2020-12-31")).values)], ts_gesla.sel(time=slice("1979-01-01","2020-12-31")).values[~np.isnan(ts_gesla.sel(time=slice("1979-01-01","2020-12-31")).values)])
            ax4.text(tmp,-tmp,f'corr={np.round(corr.statistic,3)}', verticalalignment ='bottom', horizontalalignment ='right', fontsize = 14, fontweight ='bold');
            ax4.set_xlabel('GTSM-ERA5 [m]'); ax4.set_ylabel('GESLA [m]'); ax4.set_title('1979-2020')

            figname = f'validation_{str(int(ts_gtsm.stations.values)).zfill(5)}_{ds_ges.station_name.values[ids_ges[ss]]}_detrended.png' 
            fig.savefig(f'{dir_out}/{figname}')
            plt.clf()

        del ts_gtsm, ts_gesla

    # export data

    df_stats = pd.DataFrame({
            "station": stations[ids_gtsm],
            "station_name": ds_ges.station_name.values[ids_ges],
            "bias_day_1950_1978": all_bias_day[:,0],
            "bias_day_1979_2020": all_bias_day[:,1],
            "mae_day_1950_1978": all_mae_day[:,0],
            "mae_day_1979_2020": all_mae_day[:,1],
            "mape_day_1950_1978": all_mape_day[:,0],
            "mape_day_1979_2020": all_mape_day[:,1],
            "rmse_day_1950_1978": all_rmse_day[:,0],
            "rmse_day_1979_2020": all_rmse_day[:,1],
            "corr_day_1950_1978": all_corr_day[:,0],
            "corr_day_1979_2020": all_corr_day[:,1],

            "bias_month_1950_1978": all_bias_mon[:,0],
            "bias_month_1979_2020": all_bias_mon[:,1],
            "mae_month_1950_1978": all_mae_mon[:,0],
            "mae_month_1979_2020": all_mae_mon[:,1],
            "mape_month_1950_1978": all_mape_mon[:,0],
            "mape_month_1979_2020": all_mape_mon[:,1],
            "rmse_month_1950_1978": all_rmse_mon[:,0],
            "rmse_month_1979_2020": all_rmse_mon[:,1],
            "corr_month_1950_1978": all_corr_mon[:,0],
            "corr_month_1979_2020": all_corr_mon[:,1],

            "bias_year_1950_1978": all_bias_yr[:,0],
            "bias_year_1979_2020": all_bias_yr[:,1],
            "mae_year_1950_1978": all_mae_yr[:,0],
            "mae_year_1979_2020": all_mae_yr[:,1],
            "mape_year_1950_1978": all_mape_yr[:,0],
            "mape_year_1979_2020": all_mape_yr[:,1],
            "rmse_year_1950_1978": all_rmse_yr[:,0],  
            "rmse_year_1979_2020": all_rmse_yr[:,1],
            "corr_year_1950_1978": all_corr_yr[:,0],
            "corr_year_1979_2020": all_corr_yr[:,1],
        })
    df_stats.to_csv(os.path.join(dir_out_stats,'statistics_validation_GTSM_vs_GESLA_detrended.csv'), index=False)


    df_stats = pd.read_csv(os.path.join(dir_out_stats,'statistics_validation_GTSM_vs_GESLA_detrended.csv'))
    
    # Make plots
    mpl.rcParams.update({'font.size': 18})
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    csize = 60
    cmap = mpl.colormaps['viridis'].resampled(20)
    cmap2 = mpl.colormaps['hot_r']#.resampled(20)
    cmap3 = mpl.colormaps['magma_r'].resampled(20)
    cmap4 = mpl.colormaps['seismic']#.resampled(20)


    fig = plt.figure(figsize=(20, 22))
    axs = [fig.add_subplot(4, 2, i+1, projection=crt.crs.Robinson()) for i in range(8)]
    fig.subplots_adjust(hspace=0.2, wspace=0.05)

    #row_titles = ['Bias', 'MAPE', 'Pearson correlation', 'RMSE']
    #row_y_positions = [0.88, 0.67, 0.47, 0.27]  # adjust to match row centres
    #for title, y in zip(row_titles, row_y_positions):
    #    fig.text(0.51, y, title, ha='center', va='bottom', fontsize=24, fontweight='bold')
    # plot bias for annual maxima
    # Remove lat/lon tick labels from all subplots

            
    ax = global_map(axs[0])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['bias_year_1950_1978'],
                   cmap='RdBu',transform=crt.crs.PlateCarree(),vmin=-0.5, vmax=0.5,
                   edgecolor='gray',zorder=15)
    ax.set_title('1950-1978',fontsize=22)

    ax = global_map(axs[1])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['bias_year_1979_2020'],
                   cmap='RdBu',transform=crt.crs.PlateCarree(),vmin=-0.5, vmax=0.5,
                   edgecolor='gray',zorder=15)
    ax.set_title('1979-2020',fontsize=22)
    cbar = fig.colorbar(bs, ax=[axs[0], axs[1]], orientation='horizontal', 
                        location='bottom', pad=0.1, shrink=0.5, aspect=40)
    cbar.set_label('Bias [m]',fontsize=20)

    ax = global_map(axs[2])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['mape_year_1950_1978']*100,
                   cmap=cmap2,transform=crt.crs.PlateCarree(),vmin=0, vmax=50,
                   edgecolor='gray',zorder=15)
    ax.set_title('1950-1978',fontsize=22)

    ax = global_map(axs[3])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['mape_year_1979_2020']*100,
                   cmap=cmap2,transform=crt.crs.PlateCarree(),vmin=0, vmax=50,
                   edgecolor='gray',zorder=15)
    ax.set_title('1979-2020',fontsize=22)
    cbar = fig.colorbar(bs, ax=[axs[2], axs[3]], orientation='horizontal',
                        location='bottom', pad=0.1, shrink=0.5, aspect=40)
    cbar.set_label('MAPE [%]',fontsize=20)

    # plot Pearson corr. coeff.
    ax = global_map(axs[4])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['corr_year_1950_1978'],
                   cmap=cmap3,transform=crt.crs.PlateCarree(),vmin=0.4, vmax=1,
                   edgecolor='gray',zorder=15)
    ax.set_title('1950-1978',fontsize=22)

    ax = global_map(axs[5])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['corr_year_1979_2020'],
                   cmap=cmap3,transform=crt.crs.PlateCarree(),vmin=0.4, vmax=1,
                   edgecolor='gray',zorder=15)
    ax.set_title('1979-2020',fontsize=22)
    cbar = fig.colorbar(bs, ax=[axs[4], axs[5]], orientation='horizontal', 
                        location='bottom', pad=0.1, shrink=0.5, aspect=40)
    cbar.set_label('Pearson correlation [-]',fontsize=20)

    ax = global_map(axs[6])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['rmse_year_1950_1978'],
                   cmap='viridis',transform=crt.crs.PlateCarree(),vmin=0, vmax=1,
                   edgecolor='gray',zorder=15)
    ax.set_title('1950-1978',fontsize=22)

    ax = global_map(axs[7])
    bs = ax.scatter(x=ds_gtsm.station_x_coordinate.values,
                   y=ds_gtsm.station_y_coordinate.values,
                   s=csize,c=df_stats['rmse_year_1979_2020'],
                   cmap='viridis',transform=crt.crs.PlateCarree(),vmin=0, vmax=1,
                   edgecolor='gray',zorder=15)
    ax.set_title('1979-2020',fontsize=22)
    cbar = fig.colorbar(bs, ax=[axs[6], axs[7]], orientation='horizontal', 
                        location='bottom', pad=0.1, shrink=0.5, aspect=40)
    cbar.set_label('RMSE [m]',fontsize=20)
    figname = f'annual_maxima_comparison_detrended.jpg' 
    fig.savefig(f'{dir_out_stats}/{figname}',format='jpg',dpi=300)


 #%%
# Export station coordinates + correlations to CSV
df_corr_coords = pd.DataFrame({
    "station_name": ds_ges.station_name.values[ids_ges],
    "x": ds_ges.station_x_coordinate.values[ids_ges],
    "y": ds_ges.station_y_coordinate.values[ids_ges],
    "corr_year_1950_1978": df_stats["corr_year_1950_1978"].values,
    "corr_year_1979_2020": df_stats["corr_year_1979_2020"].values,
})

out_csv = os.path.join(dir_out_stats, "gesla_station_coords_and_corr_year_detrended.csv")
df_corr_coords.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")
# %%
