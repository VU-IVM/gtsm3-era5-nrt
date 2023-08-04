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

# Author: Sanne Muis
# Contact: sanne.muis@deltares.nl
# Date created: July 2023
# Remarks: gesla

# +
#import modin.pandas as pd
import pandas as pd
from global_map import global_map
from eva import read_gesla, process_gesla, detrend, stats, compute_eva 
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
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
# -

if __name__ == "__main__":   
    settings = {'yearmin': 1950,
                'yearmax': 2022,
                'record_length': 30,
                'prc_data_missing': 25}
    prcts = [0.90,0.95,0.99,0.999]
    rps =[2,10,50,100]
    odir = 'output'
    cpu_num= 18
    # load and process GESLA timeseries  
    obj_gesla, filenames = read_gesla('data/GESLA')
    ofile = 'ds_gesla_%s_%s_allstations_%syr_max%sprt_missing_allstations.nc' % (str(settings['yearmin']),str(settings['yearmax']),str(settings['record_length']),str(settings['prc_data_missing']))
    ofile = Path(os.path.join(odir, ofile))
    print('-------->  ofile:', ofile)
    if os.path.isfile(ofile):
        ds_gesla = xr.open_dataset(ofile)
    else:
        #results = process_gesla(obj_gesla, filenames[0], settings) # serial test
        results = (Parallel(n_jobs=cpu_num)(delayed(process_gesla)(obj_gesla, ifile, settings) for ifile in filenames))
        ds_gesla = [ds for ds in results if not len(ds.dims)==0]
        ds_gesla = xr.concat(ds_gesla, dim='stations') 
        ds_gesla.to_netcdf(ofile)
    # detrend
    ds_gesla = detrend(ds_gesla)
    # compute stats
    ds_stats = stats(ds_gesla, prcts)
    # extreme value analysis
    ofile_eva = Path(str(ofile).replace('.nc', '_eva_nc'))
    if os.path.isfile(ofile_eva):
        ds_gesla_eva = xr.open_dataset(ofile_eva)
    else:
        #results = compute_eva(ds_gesla.sea_level_detrended.isel(stations=1)) # serial test
        ds_gesla_eva = pd.concat(Parallel(n_jobs=cpu_num, prefer="threads")(delayed(compute_eva)(ds_gesla.sea_level_detrended.isel(stations=istation),istation) for istation in range(0,ds_gesla.dims['stations'])))
    

if __name__ == "__main__":  
   
    # Colormap settings
    cmap = mpl.colormaps['viridis'].resampled(20)
    vrange=[0,6]
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
        bs = ax.scatter(x=ds_stats['station_x_coordinate'],y=ds_stats['station_y_coordinate'],
                        s=50,c=ds_gesla_eva[str(rps[ii])],transform=crt.crs.PlateCarree(),
                       cmap=cmap, vmin=vrange[0], vmax=vrange[1])
        cbar = ax.cax.colorbar(bs)


