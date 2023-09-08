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
from eva import detrend, stats, compute_eva
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
import multiprocessing
from joblib import Parallel, delayed

def compute_rv_eva(yr_start,yr_end):
    settings = {'yearmin': int(yr_start),
                'yearmax': int(yr_end)}

    prcts = [0.90,0.95,0.99,0.999]
    rps =[2,10,50,100]    
    cpu_num=15
    
    # location of yearly timeseries for all stations
    sys.path.append("..")
    from path_dict import path_dict
    dir_postproc = path_dict['postproc']

    dir_ds = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly','waterlevel')
    dir_ds2 = os.path.join(dir_postproc,'timeseries-GTSM-ERA5-hourly-1979-2018','waterlevel')
    dir_out = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin"]}_{settings["yearmax"]}') #output dir
    os.makedirs(dir_out,exist_ok=True)

	#determine unique stations in ds that should be removed to match 1979-2018 dataset:
    file_nc = os.path.join(dir_postproc,f'{dir_ds}/era5_reanalysis_waterlevel_1951_01.nc')
    ds = xr.open_dataset(file_nc); ds.close()
    file_nc = os.path.join(dir_postproc,f'{dir_ds2}/reanalysis_waterlevel_hourly_1979_01_v1.nc')
    ds2 = xr.open_dataset(file_nc); ds.close()
    unique_stations = np.setxor1d(ds['stations'], ds2['stations'])
    del ds, ds2 
	
    # Loop over chunks of 1000 stations
    for bb in range(1,44):

        #open and merge data for selected stations
        for year in range(settings['yearmin'],settings['yearmax']+1):
            print('merging year ', year)
            for mnth in range(1,13):
                if ((year < 1979) | (year > 2018)):                
                    file_nc = os.path.join(dir_postproc,f'{dir_ds}/era5_reanalysis_waterlevel_{year}_{mnth:02d}.nc')
                    ds = xr.open_dataset(file_nc); ds.close()
                    ds = ds.drop_sel(stations=unique_stations)                    
                else:
                    file_nc = os.path.join(dir_postproc,f'{dir_ds2}/reanalysis_waterlevel_hourly_{year}_{mnth:02d}_v1.nc')
                    ds = xr.open_dataset(file_nc); ds.close()

	            # select specific stations
                data_sel = ds.sel(stations=slice((bb-1)*1000,bb*1000-1))
	                
                if ((year == settings['yearmin']) & (mnth == 1)):
                    ds_gtsm = data_sel
                else:
                    ds_gtsm = xr.concat([ds_gtsm,data_sel],dim="time")
					
        ds_gtsm.attrs['time_coverage_end'] = str(ds_gtsm.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())
        #print(ds_gtsm)
        
        stationlist = ds_gtsm['stations'].to_series()

        # output file for TS
        ofile = 'ds_GTSM-ERA5_%s_%s_stations_%i-%i.nc' % (str(settings['yearmin']),str(settings['yearmax']),stationlist.iloc[0],stationlist.iloc[-1])
        ofile = Path(os.path.join(dir_out, ofile))
        #ds_gtsm.to_netcdf(ofile) #save full TS for selected stations

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
        
        #print(stationlist)
        ds_gtsm_eva = pd.concat(Parallel(n_jobs=cpu_num, prefer="threads")(delayed(compute_eva)(ds_gtsm.sea_level_detrended.sel(stations=istation),istation) for istation in stationlist))
        print(ds_gtsm_eva)
        ds_gtsm_eva.to_csv(ofile_eva)

        # Colormap settings
        cmap = mpl.colormaps['viridis'].resampled(20)
        vrange=[0,6]
        # Plot
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
            bs = ax.scatter(x=ds_stats['station_x_coordinate'],y=ds_stats['station_y_coordinate'],
                            s=50,c=ds_gtsm_eva[str(rps[ii])],transform=crt.crs.PlateCarree(),
                        cmap=cmap, vmin=vrange[0], vmax=vrange[1])
            cbar = ax.cax.colorbar(bs)
            figname = 'EVA_map_%s_%s_stations_%i-%i.png' % (str(settings['yearmin']),str(settings['yearmax']),stationlist.iloc[0],stationlist.iloc[-1])
        fig.savefig(f'{dir_out}/{figname}')
        print('Done!')

if __name__ == "__main__":   
    if len(os.sys.argv)>1:
        yr_start = os.sys.argv[1]
        yr_end = os.sys.argv[2]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate start year, second - end year.')
    compute_rv_eva(yr_start,yr_end)

