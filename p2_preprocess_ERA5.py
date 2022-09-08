#!/usr/bin/env python
  
def convert2FM(idir,tstart,tend):
  """
  This function converts a ERA5 file from CDS to FM format
  """
  # import modules
  import os
  import xarray as xr
  import datetime as dt
  import numpy as np
  import pandas as pd
  import glob
  # debug
  debug=False
  if debug==True:
    tstart=dt.datetime.strptime('22/12/2021', '%d/%m/%Y')
    tend=dt.datetime.strptime('01/01/2023', '%d/%m/%Y')
    idir="/gpfs/work1/0/einf3499/meteo_ERA5"
  # import dataset
  no_files = len(glob.glob(os.path.join(idir,"ERA5_CDS_atm_%s*.nc" %(tstart.year))))
  print('no. of files:', no_files)
  ds = xr.open_mfdataset(os.path.join(idir,"ERA5_CDS_atm_%s*.nc" %(tstart.year)))
  spinup = dt.datetime.strptime('12-%s' %(tstart.year-1), "%m-%Y") 
  ds_spinup = xr.open_mfdataset(os.path.join(idir,"ERA5_CDS_atm_%s*.nc" %(spinup.strftime('%Y-%m'))))
  ds_merged = xr.combine_by_coords([ds, ds_spinup],compat='equals',combine_attrs='override').sel(time=slice(tstart,tend))
  # copy latitude and create new longitude
  lats = ds['latitude'][:]
  lons = ds['longitude'][:]
  part1 = (ds.longitude>178) #move 180:360 part to -180:0 so field now runs from longitute -180 to 180
  part2 = (ds.longitude<182) #take a bit of overlap to avoid interpolation issues at edge
  lons_new=np.hstack((lons[part1]-360,lons[part2]))
  # times
  reftime = dt.datetime(1900,1,1)
  times = ds['time'][:]
  ntimes = len(times)
  dates = pd.to_datetime(pd.DataFrame({'year':times.dt.year, 'month':times.dt.month, 'day':times.dt.day, 'hour': times.dt.hour, 'minute':times.dt.minute, 'second':times.dt.second}))
  rel_times = (dates-reftime) / np.timedelta64(1, 'h')
  # copy to new dataset    
  coords =  {'latitude': lats, 'longitude': lons_new, 'time': rel_times[0]}
  # v10
  datasets = []
  for i in range(0,10):#len(rel_times)):
    v10 = ds.v10.isel(time=i).values
    var_new=np.hstack((v10[:,part1],v10[:,part2]))
    da = xr.DataArray(var_new, coords=coords, dims=['latitude', 'longitude'])
    da.name = 'v10'
    dat = xr.concat([da], 'time')
    datasets.append(dat)
  ds_v10_merged = xr.concat(datasets, dim='time') 
  print(ds_v10_merged)
  # u10
  datasets = []
  for i in range(0,10):#len(rel_times)):
    u10 = ds.u10.isel(time=i).values
    var_new=np.hstack((u10[:,part1],u10[:,part2]))
    da = xr.DataArray(var_new, coords=coords, dims=['latitude', 'longitude'])
    da.name = 'u10'
    dat = xr.concat([da], 'time')
    datasets.append(dat)
  ds_u10_merged = xr.concat(datasets, dim='time') 
  # mslp
  datasets = []
  for i in range(0,10):#len(rel_times)):
    msl = ds.msl.isel(time=i).values
    var_new=np.hstack((msl[:,part1],msl[:,part2]))
    da = xr.DataArray(var_new, coords=coords, dims=['latitude', 'longitude'])
    da.name = 'msl'
    dat = xr.concat([da], 'time')
    datasets.append(dat)
  ds_msl_merged = xr.concat(datasets, dim='time') 
  # merge and set attributes
  ds_all = xr.merge([ds_v10_merged,ds_u10_merged,ds_msl_merged]) 
  # change attributes    
  ds_all.u10.attrs['standard_name'] = 'eastward_wind'
  ds_all.u10.attrs['long_name'] = '10 metre U wind component'
  ds_all.u10.attrs['units'] = 'm s**-1'
  ds_all.v10.attrs['standard_name'] = 'northward_wind'
  ds_all.v10.attrs['long_name'] = '10 metre V wind component'
  ds_all.v10.attrs['units'] = 'm s**-1'
  ds_all.msl.attrs['standard_name'] = 'air_pressure'
  ds_all.msl.attrs['long_name'] = 'Mean sea level pressure'
  ds_all.msl.attrs['standard_name'] = 'Pa'
  # user defined scaling factor and offset
  scalef_msl=float(1)
  scalef_u10=float(0.01)
  scalef_v10=float(0.01) 
  offset_msl=float(100000)
  offset_u10=float(0)
  offset_v10=float(0)
  encoding = {}
  encoding['msl'] = {'scale_factor': scalef_msl, 'add_offset': offset_msl}
  encoding['u10'] = {'scale_factor': scalef_u10, 'add_offset': offset_u10}
  encoding['v10'] = {'scale_factor': scalef_v10, 'add_offset': offset_v10}
  ofile=os.path.join(idir.replace('meteo_ERA5','meteo_ERA5_fm'),"ERA5_CDS_atm_%s_%s.nc" %(dt.datetime.strftime(tstart, "%Y-%m-%d"),dt.datetime.strftime(tend, "%Y-%m-%d")))
  print("Saving file as:",ofile)
  ds_all.to_netcdf(ofile, encoding=encoding)

if __name__ == "__main__":
    # read input arguments
    import os
    import datetime
    if len(os.sys.argv)>0:
      tstart=os.sys.argv[2]
      tend=os.sys.argv[3]
      outdir=os.sys.argv[1]        
    else:
      raise RuntimeError('No arguments were provided\nFirst argument should indicate startdate as "%Y-%m-%d".\n Second argument for outdir. Script will download monthly files per day')
    convert2FM(tstart,tend,outdir)