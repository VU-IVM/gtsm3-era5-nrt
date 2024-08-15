#!/usr/bin/env python

import os
import calendar
import cdsapi
import sys
sys.path.append("..")
from path_dict import path_dict

def download_era5(yr, mnth):
    yr=int(yr)
    mnth=int(mnth)
    outdir = path_dict['meteo_raw']
    # find times for monthly downlods
    numdays = calendar.monthrange(yr,mnth)[1] #returns weekday of first day of month and the number of days in the month
    # daily download
    print ("######### ERA-5 from CDS #########")
    os.makedirs(outdir,exist_ok=True)
    # I/O - download the data
    for day in range(1,numdays+1):
        date_str = f'{yr}-{mnth:02d}-{day:02d}'
        print (f'getting data for {date_str}')
        targetfile = os.path.join(outdir,f"ERA5_CDS_atm_{date_str}.nc")
        if os.path.exists(targetfile):
            print('already downloaded')
            pass
        if os.path.isfile(targetfile)==False: 
            c = cdsapi.Client()
            c.retrieve('reanalysis-era5-single-levels',
                {'product_type':'reanalysis',
                'data_format':'netcdf',
                'download_format':'unarchived',
                'variable':['10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure'],
                'year':yr,
                'month':[f"{mnth:02}"],
                'day':day,
                'time':['00:00','01:00','02:00','03:00','04:00','05:00',
                        '06:00','07:00','08:00','09:00','10:00','11:00',
              	        '12:00','13:00','14:00','15:00','16:00','17:00',
                        '18:00','19:00','20:00','21:00','22:00','23:00']
                },targetfile)


if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        yr = os.sys.argv[1]
        mnth = os.sys.argv[2]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year. Second argument should indicate month. Script will download monthly files per day')
    download_era5(yr,mnth)
