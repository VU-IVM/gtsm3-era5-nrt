#!/usr/bin/env python

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

  tstart = dt.datetime.strptime(yr, "%Y").date()
  
  # define dates
  spinup_period = [1,1,15,0] # imposed 1 day zero, 1 day transition, 15 days spinup --> compute days of spinup required
  date_start_zero = dt.datetime(tstart.year,1,1)-dt.timedelta(days=int(np.sum(spinup_period[0:3]))) #eg 15 dec
  date_start_transition = date_start_zero+dt.timedelta(days=spinup_period[0]) #eg 16 dec
  date_start_spinup = date_start_zero+dt.timedelta(days=spinup_period[0]+spinup_period[1]) #eg 17 dec
  date_end = dt.datetime(tstart.year+1,1,1)

  script_tstart = dt.datetime.now()

  varkey_list = ['msl','u10','v10'] #charnock, mean_sea_level_pressure, 10m_u_component_of_neutral_wind, 10m_v_component_of_neutral_wind

  drop_variables = None
  preprocess = None
  rename_variables = None

  # Create output folder
  dir_output = input_dir.replace('meteo_ERA5','meteo_ERA5_fm')
  if not os.path.exists(dir_output):
      os.makedirs(dir_output)

  #generating file list
  dir_data = os.path.join(input_dir,f'ERA5_CDS_atm_*-*-*.nc')

   
  file_list = glob.glob(dir_data)

  #print(f'opening multifile dataset of {len(file_list)} files matching "{fn_match_pattern}" (can take a while with lots of files)')
  #open_mfdataset when dropping x/y since both varname as dimname, which is not possible in xarray
  data_xr = xr.open_mfdataset(dir_data,
                              drop_variables=drop_variables, #necessary since dims/vars with equal names are not allowed by xarray, add again later and requested matroos to adjust netcdf format.
                              parallel=True, #speeds up the process
                              preprocess=preprocess,
                              chunks={'time':1}).sel(time=slice(date_start_zero,date_end)); data_xr.close()
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
      raise Exception(f'ERROR: date_start_zero="{date_start_zero}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
  if not date_end in times_pd.index:
      raise Exception(f'ERROR: date_end="{date_end}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')

  #get actual tstart/tstop strings from times_pd
  tstart_str = times_pd.index[0].strftime("%Y%m%d")
  tstop_str = times_pd.index[-1].strftime("%Y%m%d")
  
  convert_360to180 = (data_xr['longitude'].to_numpy()>180).any()
  if convert_360to180: #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
    data_xr.coords['longitude'] = (data_xr.coords['longitude'] + 180) % 360 - 180
    data_xr = data_xr.sortby(data_xr['longitude'])
  
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
  data_xr = data_xr.where(data_xr.time >= np.datetime64(date_start_spinup), 0) # values to zero for initalization SLR
  bool = (data_xr.time > np.datetime64(date_start_transition)) & (data_xr.time < np.datetime64(date_start_spinup)) # select transition period
  data_xr = data_xr.where(~bool, drop=True) # drop times for transition period
  
      #field_zerostart = data_xr.isel(time=[0,0])*0 #two times first field, set values to 0
#      field_zerostart['time'] = [times_pd.index[0]-dt.timedelta(days=2),times_pd.index[0]-dt.timedelta(days=1)] #TODO: is one zero field not enough? (is replacing first field not also ok? (results in 1hr transition period)
#      data_xr = xr.concat([field_zerostart,data_xr],dim='time')#.sortby('time')

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