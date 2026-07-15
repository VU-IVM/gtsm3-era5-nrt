# Get surge timeseries from a selection of GESLA3 stations
#Author: n-aleksandrova

#%%
# General libraries for file paths, data extraction, etc
from glob import glob
import os

# Libraries for working with multi-dimensional arrays
import numpy as np
import xarray as xr
import pandas as pd

# Libraries for plotting and visualising data
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Libraries for dealing with time variables
from datetime import datetime

# library for tidal analysis
import hatyan

#%%
ds_ges = xr.open_dataset(r'p:\archivedprojects\11210221-gtsm-reanalysis\GESLA\ds_gesla_1950_2022_allstations_50yr_max25prt_missing.nc')

stationlist = pd.read_csv(r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\validation_stats\GESLA3_noduplicates_max_25perc_missing_both_periods.csv', index_col=0)

dir_fig = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\tidal_analysis\figures'
dir_data = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\tidal_analysis\data'

# %%
# Select only stations from the list
mask = ds_ges.station_name.isin(stationlist)
ds_ges_sel = ds_ges.isel(stations=mask.values)

# drop duplicate stations (keep first occurrence)
_, unique_idx = np.unique(ds_ges_sel.station_name.values, return_index=True)
ds_ges_sel = ds_ges_sel.isel(stations=unique_idx)

#%% 
#tidal analysis
const_list = hatyan.get_const_list_hatyan("year")

for ii, stname in enumerate(ds_ges_sel.station_name.values):

    # if ii < 60:
    #     continue
    stname = str(stname)
    print(f"{ii}: {stname}")
    st_x = ds_ges_sel.sel(station_name=stname).station_x_coordinate.values
    st_y = ds_ges_sel.sel(station_name=stname).station_y_coordinate.values
    
    ts_ges = pd.DataFrame(
                {"values": ds_ges_sel["sea_level"].sel(station_name=stname).values},
                index=pd.to_datetime(ds_ges_sel.time.values),
            )

    ts_ges = ts_ges[~ts_ges.index.duplicated()]
    ts_ges = ts_ges.dropna()

    ts_ges = hatyan.resample_timeseries(
                    ts=ts_ges,
                    timestep_min=60,
                )
    
    available_start = ts_ges.index.min()
    available_end = ts_ges.index.max()
    times_ext = [available_start, available_end]

    print(f"    Available period: {available_start} to {available_end}")

    ts_ges = hatyan.crop_timeseries(
        ts=ts_ges,
        times=slice(times_ext[0], times_ext[1]),
    )

    # Expected hours per year (accounts for leap years)
    expected = ts_ges.resample('YE').size().apply(lambda n: n)  # actual slots after resampling
    valid = ts_ges['values'].resample('YE').count()           # non-NaN count

    # Keep only years with >= 50% valid data
    good_years = valid.index.year[valid.values / expected.values >= 0.8]

    ts_ges_clean = ts_ges[ts_ges.index.year.isin(good_years)]
    
    comp_mean, comp_allperiods = hatyan.analysis(ts=ts_ges_clean, const_list=const_list, 
                                                nodalfactors=True, return_allperiods=True,
                                                fu_alltimes=True, analysis_perperiod='Y')
    
    fig,(ax1,ax2) = hatyan.plot_components(comp=comp_mean, comp_allperiods=comp_allperiods)
    ax1.legend([]); ax2.legend([])
    fig.savefig(os.path.join(dir_fig, f"{stname}_tidal_analysis.png"), dpi=200)
    plt.close()

    ts_prediction = hatyan.prediction(comp=comp_mean, times=ts_ges.index)
    fig, (ax1,ax2) = hatyan.plot_timeseries(ts=ts_prediction, ts_validation=ts_ges)
    ax1.legend(['prediction','measurement','difference','mean of prediction'])
    ax2.set_ylim(-1,1)
    fig.savefig(os.path.join(dir_fig, f"{stname}_prediction.png"), dpi=200)
    plt.close()

    # add surge to the original dataset
    ts_ges = ts_ges.rename(columns={"values": "sea_level"})
    ts_ges['tide'] = ts_prediction['values']
    ts_ges['surge'] = ts_ges['sea_level'] - ts_ges['tide']

    # Convert to xarray and save to netcdf
    ds_out = ts_ges.to_xarray().rename({'index': 'time'})
    ds_out = ds_out.assign_coords(station_name=stname, 
                                  station_x_coordinate=st_x, 
                                  station_y_coordinate=st_y)
    ds_out['time'].encoding = {'units': 'minutes since 1900-01-01', 'calendar': 'gregorian'}
    ds_out.to_netcdf(os.path.join(dir_fig, f"{stname}_tide_surge.nc"))

    del ts_ges, ts_ges_clean, comp_mean, comp_allperiods, ts_prediction, ds_out

#%% make checkplot of surge for year 2010
year = 2010
dir_gtsm = r'p:\archivedprojects\11210221-gtsm-reanalysis\GTSM-ERA5-E_dataset'

files_gtsm = glob(os.path.join(dir_gtsm, 'waterlevel', f"reanalysis_waterlevel_hourly_{year}_*_v3.nc"))

ds_wl = xr.open_mfdataset(files_gtsm)
ds_wl = ds_wl.chunk({'time': 720, 'stations': 2000})
ds_wl

files_gtsm_sur = [x.replace('waterlevel','surge') for x in files_gtsm]
ds_sur = xr.open_mfdataset(files_gtsm_sur)
ds_sur = ds_sur.chunk({'time': 720, 'stations': 2000}) 

for ii, stname in enumerate(ds_ges_sel.station_name.values):

    # if ii < 60:
    #     continue
    stname = str(stname)
    print(f"{ii}: {stname}")

    # open GESLA data that includes tide and surge
    ds_ges_st = xr.open_dataset(os.path.join(dir_fig, f"{stname}_tide_surge.nc"))
    ds_ges_st = ds_ges_st.sel(time=slice(f"{year}-01-01", f"{year}-12-31")).load()
    ds_ges_st['surge'] = ds_ges_st['surge'] - ds_ges_st['surge'].mean()
    # resample to daily max
    ds_ges_st = ds_ges_st.resample(time='1D').max()

    ges_x = ds_ges_st.station_x_coordinate.values
    ges_y = ds_ges_st.station_y_coordinate.values

    # find the index of the GTSM output point nearest to the tide gauge location
    abslat = np.abs(ds_wl.station_y_coordinate.values-ges_y)
    abslon = np.abs(ds_wl.station_x_coordinate.values-ges_x)
    c = np.maximum(abslon, abslat)
    ([iloc]) = np.where(c == np.min(c))
    station = ds_wl.stations.values[iloc[0]]

    ds_wl_st = ds_wl.sel(stations=station).load()
    ds_sur_st = ds_sur.sel(stations=station).load()
    ds_sur_st['surge'] = ds_sur_st['surge'] - ds_sur_st['surge'].mean()
    ds_sur_st = ds_sur_st.resample(time='1D').max()

    # reindex GESLA to GTSM time axis, filling missing timesteps with NaN
    ds_ges_st = ds_ges_st.reindex(time=ds_sur_st.time)

    fig,ax = plt.subplots(figsize=(10,5), nrows=2, ncols=1, sharex=True)
    ds_ges_st['surge'].plot(ax=ax[0], label='GESLA3')
    ds_sur_st['surge'].plot(ax=ax[0], alpha=0.7, label='GTSM')
    ax[0].grid(); ax[0].set_xlabel('');
    ax[0].legend(); ax[0].set_title('')
    ax[0].set_xlim([pd.to_datetime(f"{year}-01-01"), pd.to_datetime(f"{year}-12-31")])
    (ds_sur_st['surge'] - ds_ges_st['surge']).plot(ax=ax[1], label='Difference (GTSM - GESLA)', color='k')
    ax[1].grid()
    ax[1].legend(); ax[1].set_title('')
    # add text box with correlation coefficient and NRMSE
    # use only timesteps where both signals are valid (no NaNs)
    m = ds_sur_st['surge'].values
    o = ds_ges_st['surge'].values
    valid = ~(np.isnan(m) | np.isnan(o))
    corr = np.corrcoef(m[valid], o[valid])[0,1]
    nrmse = np.sqrt(np.mean((m[valid] - o[valid])**2)) / np.sqrt(np.mean(o[valid]**2))
    textstr = f'Correlation: {corr:.2f}\nNRMSE:      {nrmse:.2f}'
    ax[1].text(0.98, 0.04, textstr, transform=ax[1].transAxes, fontsize=12, 
               verticalalignment='bottom', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=1))
    fig.suptitle(f"Surge comparison for {stname} ({year})")
    fig.savefig(os.path.join(dir_fig, f"{stname}_surge_comparison_GTSM_{year}.png"), dpi=200)
    plt.close()

    ds_ges_st.close(); ds_wl_st.close(); ds_sur_st.close(); del ds_ges_st, ds_wl_st, ds_sur_st, fig, ax



    





# %%
