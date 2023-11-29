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
from pyextremes import EVA
from pyextremes import plot_threshold_stability
import pandas as pd

def detrend(ds_gesla: xr.DataArray, plot = False):
  ''' remove annual means and overall mean '''
  ds_gesla = ds_gesla.assign_coords(year=ds_gesla.time.dt.strftime("%Y"))
  ds = (ds_gesla.groupby("year") - ds_gesla.groupby("year").mean("time"))
  ds_gesla['sea_level_detrended'] = ds['sea_level'] - ds['sea_level'].mean()
  ds_gesla['sea_level_detrended'] = ds_gesla['sea_level_detrended'].chunk({"time": -1, "stations": "auto"})
  if plot == True:
      fig, axs = plt.subplots(nrows=2)
      ds_gesla.sea_level.plot.line(x='time',ax=axs[0], add_legend=False)   
      ds_gesla.sea_level_detrended.plot.line(x='time',ax=axs[1],add_legend=False)  
  return ds_gesla
    
def stats(ds_gesla: xr.Dataset ,prcts):
  ds_gesla_stats = ds_gesla.sea_level_detrended.quantile(prcts, dim=('time')) 
  return ds_gesla_stats

def compute_eva_pyextremes(var,istation,rps,dir_out):
    var = var.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
    varth = var.quantile(0.99)
    model = EVA(var)
    model.get_extremes("POT", threshold=varth, r="72H"); 
    model.fit_model(distribution='genpareto',model='MLE')

    # save a diagnostic figure
    #if (istation%10) == 0:
    #    fig,ax = model.plot_diagnostic(alpha=0.95)
    #    ax[0].set_ylabel('Total water level [m]')
    #    figname = f'EVA_station_{str(istation).zfill(5)}_diagnostic_plot.png' 
    #    fig.savefig(f'{dir_out}/checkplots/{figname}'); plt.close(fig)
    
    #    fig2,ax2 = model.plot_return_values(return_period=rps, alpha=0.95); 
    #    figname = f'EVA_station_{str(istation).zfill(5)}_RVplot.png' 
    #    ax2.set_ylabel('Total water level [m]')
    #    fig2.savefig(f'{dir_out}/checkplots/{figname}'); plt.close(fig2)
        
    summary = model.get_summary(return_period=rps,alpha=0.95)
    
    summary1 = summary.drop(columns=['lower ci','upper ci'])
    for i in range(0,len(rps)):
        summary1 = summary1.rename(index={summary1.index[i]:f'{rps[i]}_bf'})    
    summary2 = summary.drop(columns=['return value','upper ci']); summary2 = summary2.rename(columns={'lower ci':'return value'})
    for i in range(0,len(rps)):
        summary2 = summary2.rename(index={summary2.index[i]:f'{rps[i]}_lower'}) 
    summary3 = summary.drop(columns=['return value','lower ci']); summary3 = summary3.rename(columns={'upper ci':'return value'})
    for i in range(0,len(rps)):
        summary3 = summary3.rename(index={summary3.index[i]:f'{rps[i]}_higher'}) 
    summary4 = pd.DataFrame([model.distribution.mle_parameters['c'],model.distribution.mle_parameters['scale'],model.distribution.fixed_parameters['floc']],index = ['c','scale','thresh']); summary4 = summary4.rename(columns={0:'return value'})
    summary_all = pd.concat([summary1,summary2,summary3,summary4])
    
    tmp2 = pd.DataFrame([istation],index=['station'],columns=summary_all.columns)
    potmodel = pd.concat([summary_all,tmp2])    
    potmodel = potmodel.transpose().set_index('station')
    return potmodel
    
def compute_rv_eva(yr_start,yr_end,st_start,st_end,mode):
    settings = {'yearmin': int(yr_start),
                'yearmax': int(yr_end),
                'st_start': int(st_start),
                'st_end': int(st_end)}

    print(f'{yr_start}-{yr_end}, stations: {st_start}-{st_end}, {mode}')
    prcts = [0.05,0.10,0.25,0.50,0.75,0.90,0.95]
    rps = np.array([1,2,5,10,25,50,75,100])
    #cpu_num=10
    
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
    dir_out = os.path.join(dir_postproc,'EVA-GTSM-ERA5',f'period_{settings["yearmin"]}_{settings["yearmax"]}_{mode}_v2') #output dir
    os.makedirs(dir_out,exist_ok=True)
    os.makedirs(os.path.join(dir_out,'checkplots'),exist_ok=True)

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
                
    ds_gtsm = ds_gtsm.chunk({"time": -1, "stations": "auto"})
                
    ds_gtsm.attrs['time_coverage_end'] = str(ds_gtsm.time.max().dt.strftime('%Y-%m-%d %H:%M:%S').item())
    #print(ds_gtsm)
    
    stationlist = ds_gtsm['stations'].values

    # output file for TS
    ofile = 'ds_GTSM-ERA5_%s_%s_stations_%05i-%05i.nc' % (str(settings['yearmin']),str(settings['yearmax']),stationlist[0],stationlist[-1])
    ofile = Path(os.path.join(dir_out, ofile))
    #ds_gtsm.to_netcdf(ofile) #save full TS for selected stations

    # detrend
    ds_gtsm['sea_level'] = ds_gtsm['waterlevel']
    ds_gtsm = ds_gtsm.drop(['waterlevel'])
    ds_gtsm = detrend(ds_gtsm)
    ds_gtsm = ds_gtsm.drop(['sea_level'])
    ds_gtsm = ds_gtsm.chunk({"time": -1, "stations": "auto"})
    
    # compute stats
    ds_stats = stats(ds_gtsm, prcts)
    ds_stats.to_netcdf(str(ofile).replace('.nc', '_stats.nc'))
    
    # extreme value analysis
    print('Performing EVA... ')
    ofile_eva = Path(str(ofile).replace('.nc', '_eva.csv'))

    #if os.path.isfile(ofile_eva):
        #ds_gtsm_eva = xr.open_dataset(ofile_eva)
    #else:
    if not os.path.isfile(ofile_eva):
        print(ds_gtsm)
        ds_gtsm_eva = compute_eva_pyextremes(ds_gtsm.sea_level_detrended.sel(stations=stationlist[0]),stationlist[0],rps,dir_out)
        for istation in stationlist[1:]:    
            ds_gtsm_eva = pd.concat([ds_gtsm_eva,compute_eva_pyextremes(ds_gtsm.sea_level_detrended.sel(stations=istation),istation,rps,dir_out)]) # serial test
        print(ds_gtsm_eva)
        ds_gtsm_eva.to_csv(ofile_eva)
    

    #print(stationlist)
    #ds_gtsm_eva = pd.concat(Parallel(n_jobs=cpu_num, prefer="threads")(delayed(compute_eva_pyextremes)(ds_gtsm.sea_level_detrended.sel(stations=istation),istation,rps) for istation in stationlist))
    


    # Colormap settings
    #cmap = mpl.colormaps['viridis'].resampled(20)
    #vrange=[0,6]
    # Plot
    #print('Plotting... ')
    #fig = plt.figure(figsize=(26,20))
    #axes_class = (GeoAxes, dict(projection=crt.crs.Robinson()))
    #axs = AxesGrid(fig, 111, axes_class=axes_class,
    #        nrows_ncols=(2, 2),
    #        share_all=True,
    #        axes_pad=1.7,
    #       cbar_mode='each',
    #        cbar_size='3%',
    #        cbar_pad=1.0,
    #        label_mode='keep')
    #for ii in range(0,4):
    #    ax = global_map(axs[ii])
    #    bs = ax.scatter(x=ds_stats['station_x_coordinate'],y=ds_stats['station_y_coordinate'],
    #                     s=50,c=ds_gtsm_eva[str(rps[ii])],transform=crt.crs.PlateCarree(),
    #                 cmap=cmap, vmin=vrange[0], vmax=vrange[1])
    #    cbar = ax.cax.colorbar(bs)
    #    figname = 'EVA_map_%s_%s_stations_%i-%i.png' % (str(settings['yearmin']),str(settings['yearmax']),stationlist[0],stationlist[-1])
    #fig.savefig(f'{dir_out}/{figname}')
    #print('Done!')

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

