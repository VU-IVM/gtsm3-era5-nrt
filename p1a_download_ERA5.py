#!/usr/bin/env python

import os
import os.path
import click
from datetime import date, datetime
import calendar
import cdsapi

@click.command()
@click.option('--output_dir', required=True, type=str,
              help='Output directory for downloaded files',)
@click.option('--date_string', required=True, type=str,
              help='String with year and month of current month (in YYYY_MM format)',)
  
def download_era5(output_dir, date_string):
  # Get current date, for month and year information
  tdate = datetime.strptime(date_string, '%Y_%m').date()
  lastday=calendar.monthrange(tdate.year,tdate.month)[1]
  dayarray=["{:02}".format(iday) for iday in range(1,lastday+1)]     
  # Daily download
  print ("######### ERA5 from CDS  #########")
  print ('get data from ', tdate)
  print ("################################")
  for day in dayarray:
    filename = f'ERA5_CDS_atm_{tdate.year}-{tdate.month:02}-{day:02}.nc' 
    output_path = os.path.join(output_dir, filename)
    if os.path.isfile(output_path)==False: 
      c=cdsapi.Client()
      c.retrieve('reanalysis-era5-single-levels',
          {'product_type':'reanalysis',
          'format':'netcdf',
          'variable':['10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure'],
          'year':tdate.year,
          'month':["{:02}".format(tdate.month)],
          'day':day,
          'time':['00:00','01:00','02:00','03:00','04:00','05:00',
                  '06:00','07:00','08:00','09:00','10:00','11:00',
                  '12:00','13:00','14:00','15:00','16:00','17:00',
                  '18:00','19:00','20:00','21:00','22:00','23:00']
          }, output_path)
   
if __name__ == "__main__":
  download_era5()
