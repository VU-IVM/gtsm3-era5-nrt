#!/usr/bin/env python

import os
import cdsapi
import zipfile
import sys
sys.path.append("..")
from path_dict import path_dict


def download_waterlevels(yr, mnth):
    yr=int(yr)
    mnth=int(mnth)

    outdir = os.path.join(path_dict['postproc'],'timeseries-GTSM-ERA5-hourly-1979-2018')
    outdir_waterlevel = os.path.join(outdir,'waterlevel')
    outdir_surge = os.path.join(outdir,'surge')
    os.makedirs(outdir_waterlevel,exist_ok=True)
    os.makedirs(outdir_surge,exist_ok=True)

    date_str = f'{yr}-{mnth:02d}'
    
    # Monthly download
    print ("######### GTSM-tides from CDS  #########")
    print (f'getting data for {date_str}')
      
    experiment = 'reanalysis' 
      
    path_ncfile_wl = os.path.join(outdir, f'{experiment}_waterlevel_hourly_{yr}_{mnth:02}_v1.nc')
    path_ncfile_surge = os.path.join(outdir, f'{experiment}_surge_hourly_{yr}_{mnth:02}_v1.nc')

    # download water levels
    if os.path.isfile(path_ncfile_wl):
        print('already downloaded')
    else:
        path_zipfile = os.path.join(outdir_waterlevel, f'download_{date_str}.zip')
        c = cdsapi.Client()
        c.retrieve(
            'sis-water-level-change-timeseries-cmip6',
            {
                'variable': 'total_water_level',
                'experiment': experiment,
                'year': yr,
                'month':[f'{mnth:02}'],
                'temporal_aggregation': 'hourly',
            }, path_zipfile)
        # unzip and remove zipfile
        with zipfile.ZipFile(path_zipfile,"r") as zip_ref:
            zip_ref.extractall(outdir_waterlevel)
        os.remove(path_zipfile)

    # download surge
    if os.path.isfile(path_ncfile_surge):
        print('already downloaded')
    else:
        path_zipfile = os.path.join(outdir_surge, f'download_{date_str}.zip')
        c = cdsapi.Client()
        c.retrieve(
            'sis-water-level-change-timeseries-cmip6',
            {
                'format': 'zip',
                'variable': 'storm_surge_residual',
                'experiment': experiment,
                'year': yr,
                'month':[f'{mnth:02}'],
                'temporal_aggregation': 'hourly',
            }, path_zipfile)
        # unzip and remove zipfile
        with zipfile.ZipFile(path_zipfile,"r") as zip_ref:
            zip_ref.extractall(outdir_surge)
        os.remove(path_zipfile)
    
if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        yr = os.sys.argv[1]
        mnth = os.sys.argv[2]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year. Second argument should indicate month. Script will download monthly files per day')
    download_waterlevels(yr,mnth)
