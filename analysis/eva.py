import numpy as np
import pandas as pd
from scipy.stats import genpareto as gepd
idx = pd.IndexSlice
import xarray as xr
import os
from gesla import GeslaDataset
from pathlib import Path


def stats(ds_gesla: xr.Dataset ,prcts):
  ds_gesla_stats = ds_gesla.sea_level_detrended.quantile(prcts, dim=('time')) 
  return ds_gesla_stats
    
def detrend(ds_gesla: xr.DataArray, plot = False):
  ''' remove annual means and overall mean '''
  ds_gesla = ds_gesla.assign_coords(year=ds_gesla.time.dt.strftime("%Y"))
  ds = (ds_gesla.groupby("year") - ds_gesla.groupby("year").mean("time"))
  ds_gesla['sea_level_detrended'] = ds['sea_level'] - ds['sea_level'].mean()
  if plot == True:
      fig, axs = plt.subplots(nrows=2)
      ds_gesla.sea_level.plot.line(x='time',ax=axs[0], add_legend=False)   
      ds_gesla.sea_level_detrended.plot.line(x='time',ax=axs[1],add_legend=False)  
  return ds_gesla