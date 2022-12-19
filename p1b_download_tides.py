#!/usr/bin/env python

import os
import os.path
import click
from datetime import datetime
import cdsapi
import zipfile

@click.command()
@click.option('--output_dir', required=True, type=str,
              help='Output directory for downloaded files',)
@click.option('--date_string', required=True, type=str,
              help='String with year and month of current month (in YYYY_MM format)',)

def download_tides(output_dir, date_string):
  # Get current date, for month and year information
  tdate = datetime.strptime(date_string, '%Y_%m').date()   
  # Monthly download
  print ("######### GTSM-tides from CDS  #########")
  print ('get data from ', tdate)
  print ("################################")
  output_path = os.path.join(output_dir, f'future_tide_{tdate.year}_{tdate.month:02}_v1.nc')
  if os.path.isfile(output_path)==False: 
    filename = f'download_{tdate}.zip'
    output_path = os.path.join(output_dir, filename)
    c = cdsapi.Client()
    c.retrieve(
        'sis-water-level-change-timeseries-cmip6',
        {
            'format': 'zip',
            'variable': 'tidal_elevation',
            'experiment': 'future',
            'year': tdate.year,
            'month': f{tdate.month:02},
            'temporal_aggregation': '10_min',
        }, output_path)
    # unzip
    with zipfile.ZipFile(output_path,"r") as zip_ref:
      zip_ref.extractall(output_dir)
    os.remove(output_path)

if __name__ == "__main__":
    download_tides()
