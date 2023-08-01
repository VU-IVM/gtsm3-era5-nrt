# Author: Sanne Muis
# Contact: sanne.muis@deltares.nl
# Date created: July 2023
# Remarks: gesla

from gesla import GeslaDataset
import modin.pandas as pd
#import pandas as pd
from global_map import global_map
from EVA_POT_functions import peak_over_threshold, pot

import xarray as xr
import numpy as np
import os
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import cartopy as crt
import multiprocessing
from joblib import Parallel, delayed

def read_GESLA(idir):
  # source: https://gesla787883612.wordpress.com/downloads/ [Aug 2023] 
  meta_file = os.path.join(idir,"GESLA3_ALL.csv")
  data_path = os.path.join(idir,"GESLA3.0_ALL//")
  filenames = os.listdir(data_path)
  obj_gesla = GeslaDataset(meta_file=meta_file, data_path=data_path)
  return obj_gesla, filenames
  
def process_GESLA(obj_gesla, ifile, ofile, message=False)
    data, metadata = obj_gesla.file_to_pandas(filename)
    #check metadata
    if message==True:
        print('-------------------------------------------------')
        print('Processing station: ' + metadata.site_name)
        print('Overall record quality: '+metadata.overall_record_quality)
        print('Total years of observations: ', metadata.number_of_years ) 
        print('Gauge type: ' + metadata.gauge_type)      
    if metadata.overall_record_quality!= 'No obvious issues':
        continue    
    if metadata.number_of_years < record_length:
        continue
    if metadata.gauge_type != 'Coastal':
        continue  
    # filter data and resample to hourly
    data = data[data['use_flag']==1] 
    data = data.loc['{}-1-1'.format(yearmin):'{}-1-1'.format(yearmax)]
    data['time']=data.index 
    data_hourly = data.resample('H',on='time').max()
    data_hourly = data_hourly.dropna(subset='sea_level',how='all')
    data_hourly['dt'] = pd.to_datetime(data_hourly.index).strftime("%Y")
    # count missing data per year
    missing_data = data_hourly.sea_level.groupby(data_hourly.dt).transform('count')
    data_hourly = data_hourly[(missing_data/(365.25*24))> prc_data_missing]
    if plot == True:
        data_hourly.sea_level.plot()
    if len(data_hourly.sea_level) < record_length*365.25*24:
        continue  
    #  convert to xarray 
    data = data_hourly.to_xarray()
    data = data.assign_coords({"station_x_coordinate": metadata.longitude})
    data = data.assign_coords({"station_y_coordinate": metadata.latitude})
    data = data.assign_coords({"station_name": metadata.site_name})
    return data.sea_level.to_dataset())   
      
def detrend_GESLA(ds_gesla, plot = False):
  # remove annual means and overall mean
  ds_gesla = ds_gesla.assign_coords(year=ds_gesla.time.dt.strftime("%Y"))
  ds = (ds_gesla.groupby("year") - ds_gesla.groupby("year").mean("time"))
  ds_gesla['sea_level_detrended'] = ds['sea_level'] - ds['sea_level'].mean()
  if plot == True:
      fig, axs = plt.subplots(nrows=2)
      ds_gesla.sea_level.plot.line(x='time',ax=axs[0], add_legend=False)   
      ds_gesla.sea_level_detrended.plot.line(x='time',ax=axs[1],add_legend=False)  
  return ds_gesla

def stats_GESLA(ds_gesla,prcts)
  ds_gesla_stats = ds_gesla.sea_level_detrended.quantile(prcts, dim=('time')) 
  #ds_mean = ds_gesla.sea_level_detrended.mean(dim=('time'))
  #ds_std = ds_gesla.sea_level_detrended.std(dim=('time'))
  #ds_vmin = ds_gesla.sea_level_detrended.min(dim=('time'))
  #ds_vmax = ds_gesla.sea_level_detrended.min(dim=('time'))
  #ds_count = ds_gesla.sea_level_detrended.count(dim=('time'))
  return ds_gesla_stats

def compute_eva(var):
    var = var.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
    probY, potmodel = peak_over_threshold(var)
    return probY
    
if __name__ == "__main__":   
    # Settings
    yearmin = 1950
    yearmax = 2022
    record_length = 30 
    prc_data_missing = 0.75 # max 25% missing data per year
    prcts = [0.90,0.95,0.99,0.999]
    odir = 'output'
    # Parallel processing
    cpu_num = multiprocessing.cpu_count()
    print(cpu_num)
    # Load and process GESLA timeseries  
    ofile = 'ds_gesla_%s_%s_allstations_%syr_max%sprt_missing.nc' % (str(yearmin),str(yearmax),str(record_length),str((int((1-prc_data_missing)*100))))
    ofile = os.path.join(odir, ofile)
    print('--------  ofile:', ofile)
    if os.path.isfile(ofile)==True:  
      ds_gesla = xr.open_dataset(ofile)
    else:
        data = Parallel(n_jobs=cpu_num)(delayed(process_gesla(obj_gesla, ifile for ifile in filenames)
        ds_gesla = xr.concat(data, dim='stations')  
        ds_gesla.to_netcdf(ofile)
    ds_gesla = detrend_GESLA(ds_gesla)
    # Compute descriptive statistics
    ds_stats = stats_GESLA(ds_gesla, prcts)
    # Extreme value analysis
    data = Parallel(n_jobs=cpu_num)(delayed(compute_eva(ds_gesla10.sea_level_detrended.isel(stations=istation) for istation in range(ds_gesla10.dims['stations'])])
    ds_eva = xr.concat(data, dim='stations')  
    # Plot precentiles and return periods  
    if plot==True:
      # Colormap settings
      cmap = mpl.colormaps['viridis'].resampled(20)
      vrange=[0,2]
      # Plot
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
          ds = ds_prctiles.sel(quantile=prcts[ii])
          bs = ax.scatter(x=ds['station_x_coordinate'],y=ds['station_y_coordinate'],
                          s=75,c=ds,transform=crt.crs.PlateCarree(),
                         cmap=cmap, vmin=vrange[0], vmax=vrange[1])
          cbar = ax.cax.colorbar(bs)
      fig.tight_layout()
      figname = os.path.join(odir, 'figs', 'percentiles.png')
      mpl.savefig(figname)
    
