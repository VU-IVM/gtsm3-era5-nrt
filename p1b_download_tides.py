#!/usr/bin/env python

import os
import cdsapi
import zipfile
from path_dict import path_dict


def download_tides(yr, mnth):
    yr=int(yr)
    mnth=int(mnth)
    outdir = path_dict['tides_CDS']
    date_str = f'{yr}-{mnth:02d}'
    os.makedirs(outdir,exist_ok=True)
    # Monthly download
    print ("######### GTSM-tides from CDS  #########")
    print (f'getting data for {date_str}')
    experiment = 'historical' #historical/future
    path_ncfile = os.path.join(outdir, f'{experiment}_tide_{yr}_{mnth:02}_v1.nc')
    if os.path.isfile(path_ncfile):
        print('already downloaded')
    else:
        path_zipfile = os.path.join(outdir, f'download_{date_str}.zip')
        c = cdsapi.Client()
        c.retrieve(
            'sis-water-level-change-timeseries-cmip6',
            {
                'format': 'zip',
                'variable': 'tidal_elevation',
                'experiment': experiment,
                'year': yr,
                'month':[f'{mnth:02}'],
                'temporal_aggregation': '10_min',
            }, path_zipfile)
        # unzip and remove zipfile
        with zipfile.ZipFile(path_zipfile,"r") as zip_ref:
            zip_ref.extractall(outdir)
        os.remove(path_zipfile)
    
    
if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        yr = os.sys.argv[1]
        mnth = os.sys.argv[2]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year. Second argument should indicate month. Script will download monthly files per day')
    download_tides(yr,mnth)
