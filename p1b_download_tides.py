#!/usr/bin/env python

def download_tides(tstart, outdir):
  import os
  import os.path
  from datetime import datetime
  import cdsapi
  import zipfile
  
  # Get current date, for month and year information
  tdate = datetime.strptime(tstart,"%d-%m-%Y").date()   
  # Monthly download
  print ("######### GTSM-tides from CDS  #########")
  print ('get data from ', tdate)
  print ("################################")
  output_path = os.path.join(outdir, f'future_tide_{tdate.year}_{tdate.month:02}_v1.nc')
  if os.path.isfile(output_path)==False: 
    filename = f'download_{tdate}.zip'
    output_path = os.path.join(outdir, filename)
    c = cdsapi.Client()
    c.retrieve(
        'sis-water-level-change-timeseries-cmip6',
        {
            'format': 'zip',
            'variable': 'tidal_elevation',
            'experiment': 'historical',
            'year': tdate.year,
            'month':["{:02}".format(tdate.month)],
            'temporal_aggregation': '10_min',
        }, output_path)
    # unzip
    with zipfile.ZipFile(output_path,"r") as zip_ref:
      zip_ref.extractall(outdir)
    os.remove(output_path)
    
if __name__ == "__main__":
    # read input arguments
    import os
    import datetime
    if len(os.sys.argv)>0:
      tstart=os.sys.argv[1]
      outdir=os.sys.argv[2]        
    else:
      raise RuntimeError('No arguments were provided\nFirst argument should indicate startdate as "%Y-%m-%d".\n Second argument for outdir. Script will download monthly files per day')
    download_tides(tstart,outdir)