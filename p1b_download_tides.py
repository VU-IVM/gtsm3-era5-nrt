#!/usr/bin/env python

def download_tides(yr,outdir):
  # import modules
  import os
  import cdsapi
  import zipfile
  # yearly download
  print ("######### GTSM-tides from CDS  #########")
  print ('get data from ', yr)
  print ("################################")
  # I/O - download the data
  ofile = os.path.join(outdir,'download_%s.zip' % yr)
  c = cdsapi.Client()
  c.retrieve(
      'sis-water-level-change-timeseries-cmip6',
      {
          'format': 'zip',
          'variable': 'tidal_elevation',
          'experiment': 'future',
          'year': yr,
          'month': [
              '01', '02', '03',
              '04', '05', '06',
              '07', '08', '09',
              '10', '11', '12',
          ],
          'temporal_aggregation': '10_min',
      },
      ofile)
      # unzip
  with zipfile.ZipFile(ofile,"r") as zip_ref:
    zip_ref.extractall(outdir)
  os.remove(ofile)

if __name__ == "__main__":
    # read input arguments
    import os
    import datetime
    if len(os.sys.argv)>0:
      yr=os.sys.argv[1]
      outdir=os.sys.argv[2]        
    else:
      raise RuntimeError('No arguments were provided\nFirst argument should indicate year.\n Second argument for outdir. Script will download yearlyfiles per day')
    download_tides(yr,outdir)