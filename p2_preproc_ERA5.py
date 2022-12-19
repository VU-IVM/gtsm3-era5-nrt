#!/usr/bin/env python

import os
import os.path
import xarray as xr
import datetime as dt
import numpy as np
import pandas as pd
import glob

@click.command()
@click.option('--input_dir', required=True, type=str,
              help='Input directory for downloaded files',)
@click.option('--date_string', required=True, type=str,
              help='String with year and month of current month (in YYYY_MM format)',)
  
def convert2FM(input_dir,date_string):
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
  coor_dict = {
    "latitude" : {
      "standard_name" : "latitude",
      "long_name" : "latitude",
      "units" : "degrees north"},
    "longitude" : {
      "standard_name" : "longitude",
      "long_name" : "longitude",
      "units" : "degrees_east"}}
  time_unit = "hours since 1900-01-01 00:00:00" # hard coded
  # define ref/start/stop times 
  tdate = dt.datetime.strptime(date_string, '%Y_%m').date()
  spinup_period = [1,1,15,0] # imposed 1 day zero, 1 day transition, 15 days spinup --> compute days of spinup required
  date_start_zero = dt.datetime(tdate.year,tdate.month,1)-dt.timedelta(days=int(np.sum(spinup_period[0:3]))) #eg 15 dec
  date_start_transition = date_start_zero+dt.timedelta(days=spinup_period[0]) #eg 16 dec
  date_start_spinup = date_start_zero+dt.timedelta(days=spinup_period[0]+spinup_period[1]) #eg 17 dec
  date_end = dt.datetime(tdate.year,tdate.month+1,1)
  # import dataset
  input_path = os.path.join(input_dir,f'ERA5_CDS_atm_{tdate.year}-{tdate.month:02}-*') 
  ds = xr.open_mfdataset(input_path,chunks='auto',parallel=True).sel(time=slice(date_start_zero,date_end)); ds.close()
  # copy latitude and create new longitude
  lats = ds['latitude'][:]
  lons = ds['longitude'][:]
  part1 = (ds.longitude>178) #move 180:360 part to -180:0 so field now runs from longitute -180 to 180
  part2 = (ds.longitude<182) #take a bit of overlap to avoid interpolation issues at edge
  lons_new=np.hstack((lons[part1]-360,lons[part2]))
  # load data 
  for varname in var_dict.keys():
    print(varname)
    datasets = []
    for itime, time in enumerate(ds.time):
      var = ds[varname].isel(time=itime)
      var_new = np.hstack((var[:,part1],var[:,part2]))
      coords = {'latitude': lats, 'longitude': lons_new, 'time': time.values}
      da = xr.DataArray(var_new, coords=coords, dims=['latitude', 'longitude'])
      da.name = varname
      dat = xr.concat([da],'time')
      datasets.append(dat)
    ds_var_merged = xr.concat(datasets, dim='time') 
    # drop times   
    ds_var_merged = ds_var_merged.where(ds_var_merged.time >= np.datetime64(date_start_spinup), 0) # values to zero for initalization SLR
    bool = (ds_var_merged.time > np.datetime64(date_start_transition)) & (ds_var_merged.time < np.datetime64(date_start_spinup)) # select transition period
    ds_var_merged = ds_var_merged.where(~bool, drop=True) # drop times for transition period
    # set attributes + encoding
    ds_var_merged.attrs['standard_name'] = var_dict[varname]['standard_name']
    ds_var_merged.attrs['long_name'] = var_dict[varname]['long_name']
    ds_var_merged.attrs['units'] = var_dict[varname]['units']
    ds_var_merged.attrs['coordinates'] = 'longitude latitude'
    for coor in coor_dict.keys():
      ds_var_merged[coor].attrs['standard_name'] = coor_dict[coor]['standard_name']
      ds_var_merged[coor].attrs['units'] = coor_dict[coor]['units']
      ds_var_merged[coor].attrs['long_name'] = coor_dict[coor]['long_name']
    # write output
    comp = dict(scale_factor=var_dict[varname]['scale_factor'], add_offset=var_dict[varname]['offset'])  
    encoding = {varname: comp, 'time':{'units': time_unit}} 
    filename = f'ERA5_CDS_atm_{varname}_{dt.datetime.strftime(date_start_zero, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc'
    output_path = os.path.join(input_dir.replace('meteo_ERA5','meteo_ERA5_fm'), filename)
    ds_var_merged.to_netcdf(output_path, encoding=encoding, format='NETCDF4', engine='netcdf4'); ds_var_merged.close()

if __name__ == "__main__":
    convert2FM()