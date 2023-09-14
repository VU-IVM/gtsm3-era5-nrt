#!/usr/bin/env python

import os
import calendar
import cdsapi
import sys
import zipfile
sys.path.append("..")
from path_dict import path_dict

def download_gtsm(yr):
    yr=int(yr)
    outdir = path_dict['postproc']
    outdir = os.path.join(outdir,'timeseries-GTSM-ERA5-10min-1979-2018','waterlevel')

    # yearly download
    print ("######### GTSM data from CDS #########")
    os.makedirs(outdir,exist_ok=True)
    
    # I/O - download the data
    print (f'getting data for {yr}')
    targetfile = os.path.join(outdir,f"era5_reanalysis_waterlevel_{yr}.zip")
    if os.path.exists(targetfile):
        print('already downloaded')
        pass
    if os.path.isfile(targetfile)==False: 
        c = cdsapi.Client()
        c.retrieve('sis-water-level-change-timeseries-cmip6',
            {'experiment':'reanalysis',
            'format':'zip',
            'variable':'total_water_level',
            'temporal_aggregation': '10_min',
            'year':yr,
            'month':['01','02','03','04','05','06','07','08','09','10','11','12'],
            },targetfile)
        
    # unzip and remove zipfile
    with zipfile.ZipFile(targetfile,"r") as zip_ref:
        zip_ref.extractall(outdir)
    os.remove(targetfile)


if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>0:
        yr = os.sys.argv[1]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year. Script will download yearly files.')
    download_gtsm(yr)
