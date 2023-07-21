#!/usr/bin/env python

import os
import xarray as xr
import datetime as dt
import numpy as np
import glob
from path_dict import path_dict



def convert2FM(yr):
    """
    interesting datestart is 1-1-yr, 15 days of spinup is added and two timefields with zeros and daily freq to assure zerostart
    """
    input_dir = path_dict['meteo_raw']
    add_global_overlap = True #GTSM specific: extend data beyond -180 to 180 longitude
    zerostart = True #GTSM specific: extend data with 0-value fields 1 and 2 days before date_start_spinup, so starts at 15dec in the end
    
    tstart = dt.datetime.strptime(yr, "%Y")
    date_start_spinup = tstart - dt.timedelta(days=15) # 17 dec 
    date_start_zero = date_start_spinup - dt.timedelta(days=2) # 15 dec
    date_end = dt.datetime(tstart.year+1,1,1)
    time_slice = slice(date_start_spinup,date_end)
    
    varkey_list = ['msl','u10','v10'] #charnock, mean_sea_level_pressure, 10m_u_component_of_neutral_wind, 10m_v_component_of_neutral_wind
    
    # Create output folder    
    dir_output = path_dict['meteo_fm']
    os.makedirs(dir_output,exist_ok=True)
    
    #generating file list
    dir_data = os.path.join(input_dir,'ERA5_CDS_atm_*-*-*.nc')
    file_nc_list = glob.glob(dir_data)
    file_nc_list.sort()
    print(f'>> opening multifile dataset of {len(file_nc_list)} files matching "{dir_data}" (can take a while with lots of files): ',end='')
    dtstart = dt.datetime.now()
    data_xr = xr.open_mfdataset(file_nc_list,
                                #parallel=True, #TODO: speeds up the process, but often "OSError: [Errno -51] NetCDF: Unknown file format" on WCF
                                chunks={'time':1})
    data_xr.close()
    print(f'{(dt.datetime.now()-dtstart).total_seconds():.2f} sec')
        
    data_xr = data_xr.sel(time=time_slice)
    if data_xr.get_index('time').duplicated().any():
        print('dropping duplicate timesteps')
        data_xr = data_xr.sel(time=~data_xr.get_index('time').duplicated()) #drop duplicate timesteps
    
    #check if there are times selected
    if len(data_xr.time)==0:
        raise Exception(f'ERROR: no times selected, ds_text={data_xr.time[[0,-1]].to_numpy()} and time_slice={time_slice}')

    #check if there are no gaps (more than one unique timestep)
    times_pd = data_xr['time'].to_series()
    timesteps_uniq = times_pd.diff().iloc[1:].unique()
    if len(timesteps_uniq)>1:
        raise Exception(f'ERROR: gaps found in selected dataset (are there sourcefiles missing?), unique timesteps (hour): {timesteps_uniq/1e9/3600}')
    
    #check if requested times are available in selected files (in times_pd)
    if time_slice.start not in times_pd.index:
        raise Exception(f'ERROR: time_slice_start="{time_slice.start}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    if time_slice.stop not in times_pd.index:
        raise Exception(f'ERROR: time_slice_stop="{time_slice.stop}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    
    #convert 0to360 sourcedata to -180to+180
    convert_360to180 = (data_xr['longitude'].to_numpy()>180).any()
    if convert_360to180: #TODO: make more flexible for models that eg pass -180/+180 crossing (add overlap at lon edges).
        lon_newvar = (data_xr.coords['longitude'] + 180) % 360 - 180
        data_xr.coords['longitude'] = lon_newvar.assign_attrs(data_xr['longitude'].attrs) #this re-adds original attrs
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
    # doing this drops all encoding from variables, causing them to be converted into floats. Also makes sense since 0 pressure does not fit into int16 range as defined by scalefac and offset
    if zerostart: #TODO: is one zero field not enough? (is replacing first field not also ok? (results in 1hr transition period)
        field_zerostart = xr.zeros_like(data_xr.isel(time=[0,0])) #two times first field with zeros, encoding of data_vars is dropped)
        field_zerostart['time'] = [times_pd.index[0]-dt.timedelta(days=2),times_pd.index[0]-dt.timedelta(days=1)] #this drops time var encoding
        data_xr = xr.concat([field_zerostart,data_xr],dim='time',combine_attrs='no_conflicts') #combine_attrs argument prevents attrs from being dropped
    
    #write to netcdf file
    print('writing file')
    for varname in varkey_list:
        data_xr_var = data_xr[varname]
        # original long_name and units attrs are now conserved, so do not need to be enforced here
        data_xr_var.attrs['coordinates'] = 'longitude latitude'
        filename = f'ERA5_CDS_atm_{varname}_{dt.datetime.strftime(date_start_zero, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc'
        file_out = os.path.join(dir_output, filename)
        data_xr_var.to_netcdf(file_out)


if __name__ == "__main__":
    if len(os.sys.argv)>1:
        yr=os.sys.argv[1]
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year as "yyyy".')
    tstart = dt.datetime.now()
    convert2FM(yr)
    print('time passed:',dt.datetime.now()-tstart)
    
