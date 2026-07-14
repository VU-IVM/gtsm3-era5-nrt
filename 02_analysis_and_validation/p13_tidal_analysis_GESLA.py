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

    if ii < 40:
        continue
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





