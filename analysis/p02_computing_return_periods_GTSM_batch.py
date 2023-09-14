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
from eva import detrend, stats, compute_eva, peak_over_threshold
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import xarray as xr
import numpy as np
import sys
import os
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import cartopy as crt
import glob
import multiprocessing
from joblib import Parallel, delayed

def compute_rv_eva(yr_start,yr_end,st_start,st_end,mode):
    settings = {'yearmin': int(yr_start),
                'yearmax': int(yr_end),
                'st_start': int(st_start),
                'st_end': int(st_end)}

    print(f'{yr_start}-{yr_end}, stations: {st_start}-{st_end}, {mode}')
    prcts = [0.90,0.95,0.99,0.999]
    rps =[2,10,50,100]    
    cpu_num=10
    
    # location of yearly timeseries for all stations
    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']

    if mode == '1hr':
        dir_ds = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly','waterlevel')
        dir_ds2 = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly-1979-2018','waterlevel')
    elif mode == '10min':
        dir_ds = os.path.join(dir_postproc,'timeseries-GTSM-ERA5','waterlevel')
        dir_ds2 = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-10min-1979-2018','waterlevel')
    dir_out = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin"]}_{settings["yearmax"]}_{mode}') #output dir
    os.makedirs(dir_out,exist_ok=True)

	#determine unique stations in ds that should be removed to match 1979-2018 dataset:
#    file_nc = os.path.join(dir_postproc,f'{dir_ds}/reanalysis_waterlevel*_1951_01_v1.nc'); tmp = glob.glob(file_nc); file_nc = tmp[0]
#    ds = xr.open_dataset(file_nc); ds.close(); del file_nc
#    file_nc = os.path.join(dir_postproc,f'{dir_ds2}/reanalysis_waterlevel*_1987_01_v1.nc'); tmp = glob.glob(file_nc); file_nc = tmp[0]
#    ds2 = xr.open_dataset(file_nc); ds.close()
#    unique_stations = np.setxor1d(ds['stations'], ds2['stations'])
#    del ds, ds2 

    #open and merge data for selected stations
    for year in range(settings['yearmin'],settings['yearmax']+1):
        print('merging year ', year)
        for mnth in range(1,13):
            if ((year < 1980) | (year > 2018)):     # note that 1979 run is from the extended dataset with more realistic spinup           
                file_nc = os.path.join(dir_postproc,f'{dir_ds}/reanalysis_waterlevel*_{year}_{mnth:02d}*.nc')
                tmp = glob.glob(file_nc)
                file_nc = tmp[0]
                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()
                if "station_x_coordinate" in list(ds.data_vars):
                    ds = ds.set_coords(("station_x_coordinate", "station_y_coordinate"))
                #ds = ds.drop_sel(stations=unique_stations)                    
            else:
                file_nc = os.path.join(dir_postproc,f'{dir_ds2}/reanalysis_waterlevel*_{year}_{mnth:02d}_v1.nc')
                tmp = glob.glob(file_nc)
                file_nc = tmp[0]
                ds = xr.open_dataset(file_nc,chunks={'stations': 1000}); ds.close()

            # select specific stations
            data_sel = ds.sel(stations=slice(settings['st_start'],settings['st_end']))
                
            if ((year == settings['yearmin']) & (mnth == 1)):
                ds_gtsm = data_sel
            else:
                ds_gtsm = xr.concat([ds_gtsm,data_sel],dim="time")
                
    ds_gtsm.attrs['time_coverage_end'] = str(ds_gtsm.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())
    #print(ds_gtsm)
    
    stationlist = ds_gtsm['stations'].to_series()

    # output file for TS
    ofile = 'ds_GTSM-ERA5_%s_%s_stations_%05i-%05i.nc' % (str(settings['yearmin']),str(settings['yearmax']),stationlist.iloc[0],stationlist.iloc[-1])
    ofile = Path(os.path.join(dir_out, ofile))
    ds_gtsm.to_netcdf(ofile) #save full TS for selected stations

    # detrend
    ds_gtsm['sea_level'] = ds_gtsm['waterlevel']
    ds_gtsm = ds_gtsm.drop(['waterlevel'])
    ds_gtsm = detrend(ds_gtsm)
    
    # compute stats
    ds_stats = stats(ds_gtsm, prcts)
    ds_stats.to_netcdf(str(ofile).replace('.nc', '_stats.nc'))
    
    # extreme value analysis
    print('Performing EVA... ')
    ofile_eva = Path(str(ofile).replace('.nc', '_eva.csv'))

    #if os.path.isfile(ofile_eva):
        #ds_gtsm_eva = xr.open_dataset(ofile_eva)
    #else:
        #results = compute_eva(ds_gtsm.sea_level_detrended.isel(stations=1),1) # serial test
    print(ds_gtsm)

    ####
#    fig = plt.figure(figsize=(15,6))
#    plt.plot(ds_gtsm.sea_level_detrended.isel(stations=1))
    
#    var = ds_gtsm.sea_level_detrended.isel(stations=1)
#    var = var.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
#    probY, potmodel = peak_over_threshold(var)    
#    tmp2 = pd.DataFrame([1],index=['station'],columns=potmodel.columns)
#    potmodel = pd.concat([potmodel,tmp2])    
#    potmodel = potmodel.transpose().set_index('station')

#    fig = plt.figure(figsize=(10,6))
#    plt.plot(probY['duration'],probY['emp'],'o')
#    plt.plot()
#    plt.semilogx()
#    plt.xlim(left=0)
#    plt.grid()
    ####
    
    print(stationlist)
    ds_gtsm_eva = pd.concat(Parallel(n_jobs=cpu_num, prefer="threads")(delayed(compute_eva)(ds_gtsm.sea_level_detrended.sel(stations=istation),istation) for istation in stationlist))
    print(ds_gtsm_eva)
    ds_gtsm_eva.to_csv(ofile_eva)

    # Colormap settings
    # cmap = mpl.colormaps['viridis'].resampled(20)
    # vrange=[0,6]
    # # Plot
    # print('Plotting... ')
    # fig = plt.figure(figsize=(26,20))
    # axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    # axs = AxesGrid(fig, 111, axes_class=axes_class,
    #         nrows_ncols=(2, 2),
    #         share_all=True,
    #         axes_pad=1.7,
    #         cbar_location='right',
    #         cbar_mode='each',
    #         cbar_size='3%',
    #         cbar_pad=1.0,
    #         label_mode='keep')
    # for ii in range(0,4):
    #     ax = global_map(axs[ii])
    #     bs = ax.scatter(x=ds_stats['station_x_coordinate'],y=ds_stats['station_y_coordinate'],
    #                     s=50,c=ds_gtsm_eva[str(rps[ii])],transform=crt.crs.PlateCarree(),
    #                 cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    #     cbar = ax.cax.colorbar(bs)
    #     figname = 'EVA_map_%s_%s_stations_%i-%i.png' % (str(settings['yearmin']),str(settings['yearmax']),stationlist.iloc[0],stationlist.iloc[-1])
    # fig.savefig(f'{dir_out}/{figname}')
    print('Done!')

if __name__ == "__main__":   
    if len(os.sys.argv)>4:
        yr_start = os.sys.argv[1]
        yr_end = os.sys.argv[2]
        st_start = os.sys.argv[3]
        st_end = os.sys.argv[4]
        mode = os.sys.argv[5]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate start year, second - end year. Third and fourth arguments give boundaries to which station numbers to process. Fifth specifies if the processing should be done on 1-hr or 10-min timeseries.')
    compute_rv_eva(yr_start,yr_end,st_start,st_end,mode)

