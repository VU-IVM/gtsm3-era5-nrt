#%%
import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import os
from math import sqrt, cos, radians
import numpy as np
import dfm_tools as dfmt
#%%
dir_out = r'p:\1230882-emodnet_hrsm\GTSM-ERA5\validation\gesla\annualmaxima'
ds_gtsm = xr.open_dataset(r'p:\archivedprojects\11210221-gtsm-reanalysis\GTSM-ERA5-E_dataset\EVA-GTSM-ERA5\period_1950_2022_1hr_selected_stations_GESLA.nc')
ds_ges = xr.open_dataset(r'p:\archivedprojects\11210221-gtsm-reanalysis\GESLA\ds_gesla_1950_2022_allstations_50yr_max25prt_missing.nc')

#ds_gtsm = xr.open_dataset(r'c:\Users\aleksand\OneDrive - Stichting Deltares\Documents\Projects\GTSM-ERA5\paper_GTSM-ERA5-E\analysis\gesla\period_1950_2022_1hr_selected_stations_GESLA.nc')
#ds_ges = xr.open_dataset(r'c:\Users\aleksand\OneDrive - Stichting Deltares\Documents\Projects\GTSM-ERA5\paper_GTSM-ERA5-E\analysis\gesla\ds_gesla_1950_2022_allstations_50yr_max25prt_missing.nc')

ds_gtsm.station_x_coordinate.load();
ds_gtsm.station_y_coordinate.load();

plot=0

#%%
am_rmse_all=[]
am_rmse_1951_1978_all=[]
am_rmse_1979_2022_all=[]
am_rmse_timewindow_all=[]
am_rmse_timewindow_1951_1978_all=[]
am_rmse_timewindow_1979_2022_all=[]
gtsm_ids = []


for ii, station in enumerate(ds_ges.station_name.values):
    print(f'Processing station {ii+1}/{len(ds_ges.station_name.values)}: {station}')
    ds_ges_sel = ds_ges.where(ds_ges.station_name == station, drop=True)
    if len(ds_ges_sel.stations)>1:
        ds_ges_sel = ds_ges_sel.isel(stations=0) # if multiple stations with same name, just take the first one

    # find GTSM locations corresponding to GESLA stations
    abslat = np.abs(ds_ges_sel.station_y_coordinate.item() - ds_gtsm.station_y_coordinate.values)
    abslon = np.abs(ds_ges_sel.station_x_coordinate.item() - ds_gtsm.station_x_coordinate.values)
    id=np.argmin(abslon**2 + abslat**2)
    x = (radians(ds_ges_sel.station_x_coordinate.item()) - radians(ds_gtsm.station_x_coordinate.values[id])) * cos(0.5 * (radians(ds_ges_sel.station_y_coordinate.item()) + radians(ds_gtsm.station_y_coordinate.values[id])))
    y = radians(ds_ges_sel.station_y_coordinate.item()) - radians(ds_gtsm.station_y_coordinate.values[id])
    d = 6371 * sqrt(x*x + y*y)
    if d < 10: # only consider stations less than 10 km apart
        id_gtsm = id
    else:
        id_gtsm=None
        continue

    # select GTSM data for the corresponding station and align time series
    ds_gtsm_sel = ds_gtsm.isel(stations=id_gtsm)
    ds_gtsm_sel['waterlevel'] = ds_gtsm_sel['waterlevel'] - ds_gtsm_sel['waterlevel'].sel(time=slice('1986-01-01','2005-12-31')).mean()

    ds_ges_sel = ds_ges_sel.where(~np.isnan(ds_ges_sel.sea_level), drop=True)
    ds_ges_sel['sea_level'] = ds_ges_sel['sea_level'] - ds_ges_sel['sea_level'].sel(time=slice('1986-01-01','2005-12-31')).mean()

    # Calculate annual maxima for GESLA
    ges_am_val = ds_ges_sel['sea_level'].groupby("time.year").max(dim="time")
    ges_am_time = ds_ges_sel['sea_level'].groupby("time.year").map(lambda x: x.time.isel(time=x.argmax("time")))
    annual_max_ges = xr.Dataset({
        "ges_am_value": ges_am_val,
        "ges_am_time": ges_am_time,
    })

    # Calculate annual maxima for GTSM
    gtsm_am_val = ds_gtsm_sel['waterlevel'].groupby("time.year").max(dim="time")
    gtsm_am_time = ds_gtsm_sel['waterlevel'].groupby("time.year").map(lambda x: x.time.isel(time=x.argmax("time")))
    annual_max_gtsm = xr.Dataset({
        "gtsm_am_value": gtsm_am_val,
        "gtsm_am_time": gtsm_am_time,
    })

    # Calculate annual maxima matched by time window (GESLA annual-max timestamp +/- 1.5 days)
    dt = np.timedelta64(36, "h")
    years = annual_max_ges.year

    gtsm_window_max_val = []
    gtsm_window_max_time = []

    for y in years.values:
        t0 = annual_max_ges.ges_am_time.sel(year=y).values  # GESLA annual-max timestamp for this year
        if type(t0)==np.ndarray:
            t0=t0[0]

        # GTSM values in [t0-1.5 days, t0+1.5 days]
        w = ds_gtsm_sel.waterlevel.sel(time=slice(t0 - dt, t0 + dt))

        if w.time.size == 0:
            # no GTSM data in window
            gtsm_window_max_val.append(xr.DataArray(np.nan))
            gtsm_window_max_time.append(xr.DataArray(np.datetime64("NaT")))
        else:
            i = w.argmax("time")
            gtsm_window_max_val.append(w.isel(time=i))
            gtsm_window_max_time.append(w.time.isel(time=i))

    gtsm_am_value_sametime = xr.concat(gtsm_window_max_val, dim="year").assign_coords(year=years)
    gtsm_am_time_sametime = xr.concat(gtsm_window_max_time, dim="year").assign_coords(year=years)

    # Calculate RMSE between GESLA and GTSM annual maxima
    am_rmse = sqrt(((annual_max_gtsm.gtsm_am_value - annual_max_ges.ges_am_value)**2).mean())
    am_rmse_1951_1978 = sqrt(((annual_max_gtsm.gtsm_am_value.sel(year=slice(1951,1978)) - annual_max_ges.ges_am_value.sel(year=slice(1951,1978)))**2).mean())
    am_rmse_1979_2022 = sqrt(((annual_max_gtsm.gtsm_am_value.sel(year=slice(1979,2022)) - annual_max_ges.ges_am_value.sel(year=slice(1979,2022)))**2).mean())
    #print(f"Annual Maxima RMSE: {am_rmse:.2f} m")

    am_rmse_timewindow = sqrt(((gtsm_am_value_sametime - annual_max_ges.ges_am_value)**2).mean())
    am_rmse_timewindow_1951_1978 = sqrt(((gtsm_am_value_sametime.sel(year=slice(1951,1978)) - annual_max_ges.ges_am_value.sel(year=slice(1951,1978)))**2).mean())
    am_rmse_timewindow_1979_2022 = sqrt(((gtsm_am_value_sametime.sel(year=slice(1979,2022)) - annual_max_ges.ges_am_value.sel(year=slice(1979,2022)))**2).mean())
    #print(f"Annual Maxima RMSE (matched by time window): {am_rmse_timewindow:.2f} m")

    gtsm_ids.append(id_gtsm)
    am_rmse_all.append(am_rmse)
    am_rmse_1951_1978_all.append(am_rmse_1951_1978)
    am_rmse_1979_2022_all.append(am_rmse_1979_2022)
    am_rmse_timewindow_all.append(am_rmse_timewindow)       
    am_rmse_timewindow_1951_1978_all.append(am_rmse_timewindow_1951_1978)
    am_rmse_timewindow_1979_2022_all.append(am_rmse_timewindow_1979_2022)

    # Plot location of GESLA and GTSM stations
    # fig,ax=plt.subplots(figsize=(14,8))
    # bs0 = ax.scatter(x=ds_ges_sel.station_x_coordinate.values,y=ds_ges_sel.station_y_coordinate.values,s=30,color='blue', label='GESLA'); 
    # bs = ax.scatter(x=ds_gtsm_sel.station_x_coordinate.values,y=ds_gtsm_sel.station_y_coordinate.values,s=30,color='red', label='GTSM'); 
    # ax.set_xlim([ds_ges_sel.station_x_coordinate.values-0.3, ds_ges_sel.station_x_coordinate.values+0.3])
    # ax.set_ylim([ds_ges_sel.station_y_coordinate.values-0.3, ds_ges_sel.station_y_coordinate.values+0.3])
    # dfmt.plot_coastlines(ax=ax, linewidth=0.5)
    # ax.legend()
    # ax.title.set_text(f"Location {ds_ges_sel.station_name.item()} (blue: GESLA, red: GTSM)")
    # fig.savefig(os.path.join(dir_out, 'figures', f"am_comparison_{ds_ges_sel.station_name.item()}_location.png"), dpi=300)
    # plt.close(fig)
    
    # Plot time series of GESLA and GTSM with annual maxima
    if plot:
        fig, axs = plt.subplots(figsize=(16,8), nrows=2, sharex=True, sharey=True)
        ds_ges_sel.sel(time=slice('1951-01-01','2024-12-31')).sea_level.plot(ax=axs[0], color='gray')
        axs[0].scatter(annual_max_ges.ges_am_time, annual_max_ges.ges_am_value, color='r', label='Annual Maxima')
        axs[0].set_title('GESLA')
        axs[0].grid()

        ds_gtsm_sel.sel(time=slice('1951-01-01','2024-12-31')).waterlevel.plot(ax=axs[1], color='gray')
        axs[1].scatter(annual_max_gtsm.gtsm_am_time, annual_max_gtsm.gtsm_am_value, color='b', label='Annual Maxima')
        axs[1].scatter(gtsm_am_time_sametime, gtsm_am_value_sametime, color='r', label='Annual Maxima (matched by time window)')
        axs[1].set_title('GTSM')
        axs[1].grid()

        axs[0].text(0.02, 0.2, f"RMSE: {am_rmse:.2f} m [{am_rmse_1951_1978:.2f} m, {am_rmse_1979_2022:.2f} m]\nRMSE (time window): {am_rmse_timewindow:.2f} m [{am_rmse_timewindow_1951_1978:.2f} m, {am_rmse_timewindow_1979_2022:.2f} m]", transform=axs[0].transAxes, fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        fig.suptitle(f"Annual Maxima Comparison at {ds_ges_sel.station_name.item()}", fontsize=16)
        fig.savefig(os.path.join(dir_out, 'figures', f"am_comparison_{ds_ges_sel.station_name.item()}.png"), dpi=300)

        axs[1].set_xlim([np.datetime64('1951-01-01'), np.datetime64('1978-12-31')])

        fig.savefig(os.path.join(dir_out, 'figures', f"am_comparison_{ds_ges_sel.station_name.item()}_1951-1978.png"), dpi=300)
        plt.close(fig)

    del ds_gtsm_sel, ds_ges_sel, annual_max_ges, annual_max_gtsm, gtsm_am_value_sametime, gtsm_am_time_sametime


#%%
am_rmse_1951_1978_all = np.array(am_rmse_1951_1978_all)
am_rmse_1979_2022_all = np.array(am_rmse_1979_2022_all)
am_rmse_timewindow_1951_1978_all = np.array(am_rmse_timewindow_1951_1978_all)
am_rmse_timewindow_1979_2022_all = np.array(am_rmse_timewindow_1979_2022_all)

# %%
fig, ax = plt.subplots(figsize=(14,7))
p=ax.scatter(ds_gtsm.station_x_coordinate.values[gtsm_ids], 
           ds_gtsm.station_y_coordinate.values[gtsm_ids], 
           c=am_rmse_all, cmap='viridis', s=50, edgecolor='k',
           vmin=0, vmax=0.5)
plt.colorbar(p, ax=ax, label='RMSE of Annual Maxima (m)', orientation='horizontal')
ax.set_xlabel('Longitude')  
ax.set_ylabel('Latitude')
dfmt.plot_coastlines(ax=ax, linewidth=0.5, res='l', zorder=0)

# %%
fig, ax = plt.subplots(figsize=(14,7))
p=ax.scatter(ds_gtsm.station_x_coordinate.values[gtsm_ids], 
           ds_gtsm.station_y_coordinate.values[gtsm_ids], 
           c=am_rmse_1951_1978_all-am_rmse_1979_2022_all, cmap='RdBu_r', s=50, edgecolor='k',
           vmin=-0.1, vmax=0.1)
plt.colorbar(p, ax=ax, label='RMSE of Annual Maxima difference (m)', orientation='horizontal')
ax.set_xlabel('Longitude')  
ax.set_ylabel('Latitude')
dfmt.plot_coastlines(ax=ax, linewidth=0.5, res='l',zorder=0)
# %%
fig, ax = plt.subplots(figsize=(14,7))
p=ax.scatter(ds_gtsm.station_x_coordinate.values[gtsm_ids], 
           ds_gtsm.station_y_coordinate.values[gtsm_ids], 
           c=am_rmse_timewindow_1951_1978_all-am_rmse_timewindow_1979_2022_all, cmap='RdBu_r', s=50, edgecolor='k',
           vmin=-0.1, vmax=0.1)
plt.colorbar(p, ax=ax, label='RMSE of Annual Maxima difference (m)', orientation='horizontal')
ax.set_xlabel('Longitude')  
ax.set_ylabel('Latitude')
dfmt.plot_coastlines(ax=ax, linewidth=0.5, res='l',zorder=0)
# %%
