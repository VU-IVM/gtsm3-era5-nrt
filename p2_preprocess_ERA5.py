#!/usr/bin/env python
  
def convert2FM_yearly_runs(tstart,idir):
  """
  This function converts a ERA5 file from CDS to FM format
  1. METEO-FORCING (yearly):
	- Modify for correct longitude range [-180 to 180] and overlap
	- Right var names and attributes
	- Include spinup in yearly files (10 days)
	- Set initial timesteps to zero to allow for SLR correction
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
    tstart=dt.datetime.strptime('01-01-2021', '%d-%m-%Y') 
    idir="/gpfs/work1/0/einf3499/meteo_ERA5"
  # create variables 'dictionary' for filename pattern, variable name in ncin, variable name in ncout  
  var_dict = {
    "u10" : {
      "standard_name" : "eastward_wind",
      "long_name" : "10 metre U wind component",
      "units" : "m s**-1",
      "scale_factor" : float(0.01),
      "offset" : float(0)},
    "v10" : {
      "standard_name" : "northward_wind",
      "long_name" : "10 metre V wind component",
      "units" : "m s**-1",
      "scale_factor" : float(0.01),
      "offset" : float(0)},
    "msl" : {
      "standard_name" : "air_pressure",
      "long_name" : "Mean sea level pressure",
      "units" : "Pa",
      "scale_factor" : float(1),
      "offset" : float(100000)}}              
  # define ref/start/stop times 
  spinup_period = [1,1,15,0] # imposed 1 day zero, 1 day transition, 15 days spinup
  reftime=dt.datetime(1900,1,1) #imposed
  date_start_zero = dt.datetime(tstart.year,1,1)-dt.timedelta(days=int(np.sum(spinup_period[0:3]))) #eg 15 dec
  date_start_transition = date_start_zero+dt.timedelta(days=spinup_period[0]) #eg 16 dec
  date_start_spinup = date_start_zero+dt.timedelta(days=spinup_period[0]+spinup_period[1]) #eg 17 dec
  date_end = dt.datetime(tstart.year+1,1,1)
  # import dataset
  ds = xr.open_mfdataset(os.path.join(idir,"ERA5_CDS_atm_%s*.nc" %(tstart.year)))
  ds_spinup = xr.open_mfdataset(os.path.join(idir,"ERA5_CDS_atm_%s*.nc" %(date_start_zero.strftime('%Y-%m'))))
  ds_merged = xr.combine_by_coords([ds, ds_spinup],compat='equals',combine_attrs='override').sel(time=slice(date_start_zero,date_end))
  # copy latitude and create new longitude
  lats = ds['latitude'][:]
  lons = ds['longitude'][:]
  part1 = (ds.longitude>178) #move 180:360 part to -180:0 so field now runs from longitute -180 to 180
  part2 = (ds.longitude<182) #take a bit of overlap to avoid interpolation issues at edge
  lons_new=np.hstack((lons[part1]-360,lons[part2]))
  # read times and create time
  times = ds_merged['time'][:]
  dates = pd.to_datetime(pd.DataFrame({'year':times.dt.year, 'month':times.dt.month, 'day':times.dt.day, 'hour': times.dt.hour, 'minute':times.dt.minute, 'second':times.dt.second}))
  rel_times = (dates-reftime) / np.timedelta64(1, 'h') # relative time in hours since reftime
  # copy v10, u10, msl to new dataset    
  for varname in var_dict.keys():
    datasets = []
    for itime, time in enumerate(times):
      var = ds_merged[varname].sel(time=time)
      var_new = np.hstack((var[:,part1],var[:,part2]))
      coords = {'latitude': lats, 'longitude': lons_new, 'time': time.values}
      da = xr.DataArray(var_new, coords=coords, dims=['latitude', 'longitude'])
      da.name = varname
      dat = xr.concat([da],'time')
      datasets.append(dat)
    ds_var_merged = xr.concat(datasets, dim='time') 
    ds_var_merged = ds_var_merged.where(ds_var_merged.time >= np.datetime64(date_start_spinup), 0) # values to zero for initalization SLR
    bool = (ds_var_merged.time > np.datetime64(date_start_transition)) & (ds_var_merged.time < np.datetime64(date_start_spinup)) # select transition period
    ds_var_merged = ds_var_merged.where(~bool, drop=True) # drop times for transition period
    # set attributes + encoding
    ds_var_merged.attrs['standard_name'] = var_dict[varname]['standard_name']
    ds_var_merged.attrs['long_name'] = var_dict[varname]['long_name']
    ds_var_merged.attrs['units'] = var_dict[varname]['units']
    ds_var_merged['time'].attrs['units'] = "hours since %s" % (reftime)
    comp = dict(scale_factor=var_dict[varname]['scale_factor'], add_offset=var_dict[varname]['offset'])
    encoding = {varname: comp}
    ofile=os.path.join(idir.replace('meteo_ERA5','meteo_ERA5_fm'),"ERA5_CDS_atm_%s_%s_%s.nc" %(varname, dt.datetime.strftime(tstart, "%Y-%m-%d"),dt.datetime.strftime(tend, "%Y-%m-%d")))
    print("Saving file as:",ofile)
    ds_var_merged.to_netcdf(ofile, encoding=encoding)

if __name__ == "__main__":
    # read input arguments
    import os
    import datetime as dt
    if len(os.sys.argv)>0:
      tstart=dt.datetime.strptime(os.sys.argv[1], '%Y%m%d')
      outdir=os.sys.argv[2]   
      print(tstart)     
    else:
      raise RuntimeError('No arguments were provided\nFirst argument should indicate startdate as "%Y-%m-%d".\n Second argument for outdir. Script will download monthly files per day')
    convert2FM_yearly_runs(tstart,outdir)