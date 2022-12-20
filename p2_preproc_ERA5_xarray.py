#!/usr/bin/env python

  """
  Based on script by @veenstra from Oct 2022
  Modified by @n-aleksandrova, Dec 2022
  """

#import matplotlib.pyplot as plt
#plt.close('all')


def convert2FM(yr,input_dir):

  import os
  import os.path
  import xarray as xr
  import datetime as dt
  import numpy as np
  import pandas as pd
  import glob


  add_global_overlap = True #GTSM specific: extend data beyond -180 to 180 longitude
  zerostart = True #GTSM specific: extend data with 0-value fields 1 and 2 days before all_tstart

  # define dates
  spinup_period = [1,1,15,0] # imposed 1 day zero, 1 day transition, 15 days spinup --> compute days of spinup required
  date_start_zero = dt.datetime(yr,1,1)-dt.timedelta(days=int(np.sum(spinup_period[0:3]))) #eg 15 dec
  date_start_transition = date_start_zero+dt.timedelta(days=spinup_period[0]) #eg 16 dec
  date_start_spinup = date_start_zero+dt.timedelta(days=spinup_period[0]+spinup_period[1]) #eg 17 dec
  date_end = dt.datetime(yr+1,1,1)

  script_tstart = dt.datetime.now()

  varkey_list = ['mslp','u10n','v10n'] #charnock, mean_sea_level_pressure, 10m_u_component_of_neutral_wind, 10m_v_component_of_neutral_wind

  drop_variables = None
  preprocess = None
  rename_variables = None

  # Create output folder
  dir_output = input_dir.replace('meteo_ERA5','meteo_ERA5_fm')
  if not os.path.exists(dir_output):
      os.makedirs(dir_output)

  #generating file list
  dir_data = os.path.join(input_dir,f'ERA5_CDS_atm_{yr}-*-*.nc') 
  file_list = glob.glob(dir_data)

  print(f'opening multifile dataset of {len(file_list)} files matching "{fn_match_pattern}" (can take a while with lots of files)')
  #open_mfdataset when dropping x/y since both varname as dimname, which is not possible in xarray
  data_xr = xr.open_mfdataset(dir_data,
                              drop_variables=drop_variables, #necessary since dims/vars with equal names are not allowed by xarray, add again later and requested matroos to adjust netcdf format.
                              parallel=True, #speeds up the process
                              preprocess=preprocess,
                              chunks={'time':1}).sel(time=slice(date_start_zero,date_end)); ds.close()
  print('...done')

  #rename variables
  data_xr = data_xr.rename(rename_variables)
  varkeys = data_xr.variables.mapping.keys()

  if data_xr.get_index('time').duplicated().any():
      print('dropping duplicate timesteps')
      data_xr = data_xr.sel(time=~data_xr.get_index('time').duplicated()) #drop duplicate timesteps
  times_pd = data_xr['time'].to_series()

  #check if there are times selected
  if len(times_pd)==0:
      raise Exception('ERROR: no times selected, check tstart/tstop and file_nc')

  #check if there are no gaps (more than one timestep)
  timesteps_uniq = times_pd.diff().iloc[1:].unique()
  if len(timesteps_uniq)>1:
      raise Exception(f'ERROR: gaps found in selected dataset (are there sourcefiles missing?), unique timesteps (hour): {timesteps_uniq/1e9/3600}')

  #check if requested times are available in selected files (in times_pd)
  if not date_start_zero in times_pd.index:
      raise Exception(f'ERROR: all_tstart="{all_tstart}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
  if not date_end in times_pd.index:
      raise Exception(f'ERROR: all_tstop="{all_tstop}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')

  #get actual tstart/tstop strings from times_pd
  tstart_str = times_pd.index[0].strftime("%Y%m%d")
  tstop_str = times_pd.index[-1].strftime("%Y%m%d")

  #TODO: check conversion implementation with hydro_tools\ERA5\ERA52DFM.py
  # def get_unit(data_xr_var):
  #     if 'units' in data_xr_var.attrs.keys():
  #         unit = data_xr_var.attrs["units"]
  #     else:
  #         unit = '-'
  #     return unit

  #convert Kelvin to Celcius
  # for varkey_sel in ['air_temperature','dew_point_temperature','d2m','t2m']: # 2 meter dewpoint temparature / 2 meter temperature
  #     if varkey_sel in varkeys:
  #         current_unit = get_unit(data_xr[varkey_sel])
  #         new_unit = 'C'
  #         print(f'converting {varkey_sel} unit from Kelvin to Celcius: [{current_unit}] to [{new_unit}]')
  #         data_xr[varkey_sel].attrs['units'] = new_unit
  #         data_xr[varkey_sel] = data_xr[varkey_sel] - 273.15
  #convert fraction to percentage
  # for varkey_sel in ['cloud_area_fraction','tcc']: #total cloud cover
  #     if varkey_sel in varkeys:
  #         current_unit = get_unit(data_xr[varkey_sel])
  #         new_unit = '%' #unit is soms al %
  #         print(f'converting {varkey_sel} unit from fraction to percentage: [{current_unit}] to [{new_unit}]')
  #         data_xr[varkey_sel].attrs['units'] = new_unit
  #         data_xr[varkey_sel] = data_xr[varkey_sel] * 100
  #convert kg/m2/s to mm/day
  # for varkey_sel in ['mer','mtpr']: #mean evaporation rate / mean total precipitation rate
  #     if varkey_sel in varkeys:
  #         current_unit = get_unit(data_xr[varkey_sel])
  #         new_unit = 'mm/day'
  #         print(f'converting {varkey_sel} unit from kg/m2/s to mm/day: [{current_unit}] to [{new_unit}]')
  #         data_xr[varkey_sel].attrs['units'] = new_unit
  #         data_xr[varkey_sel] = data_xr[varkey_sel] * 86400 # kg/m2/s to mm/day (assuming rho_water=1000)
  # #convert J/m2 to W/m2
  # for varkey_sel in ['ssr','strd']: #solar influx (surface_net_solar_radiation) / surface_thermal_radiation_downwards #TODO: 
  #     if varkey_sel in varkeys:
  #         current_unit = get_unit(data_xr[varkey_sel])
  #         new_unit = 'W m**-2'
  #         print(f'converting {varkey_sel} unit from J/m2 to W/m2: [{current_unit}] to [{new_unit}]')
  #         data_xr[varkey_sel].attrs['units'] = new_unit
  #         data_xr[varkey_sel] = data_xr[varkey_sel] / 3600 # 3600s/h #TODO: 1W = 1J/s, so does not make sense?
  # #solar influx increase for beta=6%
  # if 'ssr' in varkeys:
  #     print('ssr (solar influx) increase for beta=6%')
  #     data_xr['ssr'] = data_xr['ssr'] *.94

  #GTSM specific addition for longitude overlap
  if add_global_overlap: # assumes -180 to ~+179.75 (full global extent, but no overlap). Does not seem to mess up results for local models.
      if len(data_xr.longitude.values) != len(np.unique(data_xr.longitude.values%360)):
          raise Exception(f'add_global_overlap=True, but there are already overlapping longitude values: {data_xr.longitude}')
      overlap_ltor = data_xr.sel(longitude=data_xr.longitude<=-179)
      overlap_ltor['longitude'] = overlap_ltor['longitude'] + 360
      overlap_rtol = data_xr.sel(longitude=data_xr.longitude>=179)
      overlap_rtol['longitude'] = overlap_rtol['longitude'] - 360
      data_xr = xr.concat([data_xr,overlap_ltor,overlap_rtol],dim='longitude').sortby('longitude')

  #GTSM specific addition for zerovalues during spinup
  if zerostart:
      field_zerostart = data_xr.isel(time=[0,0])*0 #two times first field, set values to 0
      field_zerostart['time'] = [times_pd.index[0]-dt.timedelta(days=2),times_pd.index[0]-dt.timedelta(days=1)] #TODO: is one zero field not enough? (is replacing first field not also ok? (results in 1hr transition period)
      data_xr = xr.concat([field_zerostart,data_xr],dim='time')#.sortby('time')

  encoding = {}
  #encoding['time'] = {'units': 'hours since 1900-01-01 00:00:00'} #TODO: maybe add different reftime?
  #for varkey in list(data_xr.data_vars.keys()):
  #    encoding[varkey] = {'scale_factor':0.01,'add_offset':0} #TODO: maybe add, but not necessary since xarray uses encoding from first file and that is already quite efficient.

  #write to netcdf file
  
  print('writing file')
  
  for varname in varkey_list:
    data_xr_var = data_xr[varname]
    filename = f'ERA5_CDS_atm_{varname}_{dt.datetime.strftime(date_start_zero, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc'
    file_out = os.path.join(dir_output, filename)
    data_xr_var.to_netcdf(file_out, encoding=encoding)


#  print('loading outputfile')
#  with xr.open_dataset(file_out) as data_xr_check:
#     for varkey in data_xr_check.data_vars:
#         print(f'plotting {varkey}')
#         #continue #uncomment to skip plotting
#          fig,ax1 = plt.subplots()
#          if 'HIRLAM' in mode:
#              data_xr_check[varkey].isel(time=0).plot(ax=ax1,x='longitude',y='latitude') #x/y are necessary since coords are not 1D and dims
#          elif 'depth' in data_xr[varkey].coords:
#              data_xr_check[varkey].isel(time=0).sel(depth=0).plot(ax=ax1)
#          elif 'expver' in data_xr[varkey].coords:
#              data_xr_check[varkey].isel(time=0).isel(expver=0).plot(ax=ax1)
#          else:
#              data_xr_check[varkey].isel(time=0).plot(ax=ax1)
#          fig.savefig(file_out.replace('.nc',f'_{varkey}'))
#
#  script_telapsed = (dt.datetime.now()-script_tstart)
#  print(f'elapsed time: {script_telapsed}')

if __name__ == "__main__":
  import os
  if len(os.sys.argv)>0:
    yr=os.sys.argv[1]
    input_dir=os.sys.argv[2]        
  else:
    raise RuntimeError('No arguments were provided')
  convert2FM(yr,input_dir)