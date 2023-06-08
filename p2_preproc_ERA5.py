#!/usr/bin/env python

import os
import xarray as xr
import datetime as dt
import numpy as np
import glob

def convert2FM(yr,input_dir):
    """
    interesting datestart is 1-1-yr, 15 days of spinup is added and two timefields with zeros and daily freq to assure zerostart
    """
    add_global_overlap = True #GTSM specific: extend data beyond -180 to 180 longitude
    zerostart = True #GTSM specific: extend data with 0-value fields 1 and 2 days before date_start_spinup, so starts at 15dec in the end
    
    tstart = dt.datetime.strptime(yr, "%Y")
    date_start_spinup = tstart - dt.timedelta(days=15) # 17 dec 
    date_start_zero = date_start_spinup - dt.timedelta(days=2) # 15 dec
    date_end = dt.datetime(tstart.year+1,1,1)
    
    varkey_list = ['msl','u10','v10'] #charnock, mean_sea_level_pressure, 10m_u_component_of_neutral_wind, 10m_v_component_of_neutral_wind
    
    # Create output folder
    dir_output = input_dir.replace('meteo_ERA5','meteo_ERA5_fm')
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    
    #generating file list
    dir_data = os.path.join(input_dir,'ERA5_CDS_atm_*-*-*.nc')
    file_list = glob.glob(dir_data)
    file_list.sort()
    print(f'opening multifile dataset of {len(file_list)} files matching "{dir_data}" (can take a while with lots of files)')
    data_xr = xr.open_mfdataset(file_list,
                                #parallel=True, #TODO: speeds up the process, but often "OSError: [Errno -51] NetCDF: Unknown file format" on WCF
                                chunks={'time':1}).sel(time=slice(date_start_spinup,date_end))
    data_xr.close()
    print('...done')
    
    if data_xr.get_index('time').duplicated().any():
        print('dropping duplicate timesteps')
        data_xr = data_xr.sel(time=~data_xr.get_index('time').duplicated()) #drop duplicate timesteps
    times_pd = data_xr['time'].to_series()
    
    #check if there are times selected
    if len(times_pd)==0:
        raise Exception('ERROR: no times selected, check tstart/tstop and file_nc')
    
    #check if there are no gaps (more than one unique timestep)
    timesteps_uniq = times_pd.diff().iloc[1:].unique()
    if len(timesteps_uniq)>1:
        raise Exception(f'ERROR: gaps found in selected dataset (are there sourcefiles missing?), unique timesteps (hour): {timesteps_uniq/1e9/3600}')
    
    #check if requested times are available in selected files (in times_pd)
    if not date_start_spinup in times_pd.index:
        raise Exception(f'ERROR: date_start_spinup="{date_start_spinup}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    if not date_end in times_pd.index: #TODO: revert this comment
        raise Exception(f'ERROR: date_end="{date_end}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    
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
    #TODO: doing this drops all encoding from variables, causing them to be converted into floats. Also makes sense since 0 pressure does not fit into int16 range as defined by scalefac and offset
    #'scale_factor': 0.17408786412952254, 'add_offset': 99637.53795606793
    #99637.53795606793 - 0.17408786412952254*32768
    #99637.53795606793 + 0.17408786412952254*32767
    if zerostart:
        field_zerostart = data_xr.isel(time=[0,0])*0 #two times first field, set values to 0
        field_zerostart['time'] = [times_pd.index[0]-dt.timedelta(days=2),times_pd.index[0]-dt.timedelta(days=1)] #TODO: is one zero field not enough? (is replacing first field not also ok? (results in 1hr transition period)
        data_xr = xr.concat([field_zerostart,data_xr],dim='time',combine_attrs='no_conflicts') #combine_attrs argument prevents attrs from being dropped
    
    # create variables 'dictionary' with attrs
    var_dict = {
      "u10" : {
        "standard_name" : "eastward_wind", #TODO: missing from ERA5 dataset
        "scale_factor" : float(0.01), #TODO: scale_factor and add_offset (not offset) are not written to file, should also be written towards encoding instead of attrs. Also dtype and _FillValue encoding should be (over)written
        "offset" : float(0)},
      "v10" : {
        "standard_name" : "northward_wind",#TODO: missing from ERA5 dataset
        "scale_factor" : float(0.01),
        "offset" : float(0)},
      "msl" : {
        "standard_name" : "air_pressure", #TODO: is air_pressure_at_mean_sea_level in ERA5 dataset
        "scale_factor" : float(1),
        "offset" : float(100000)}}
    #write to netcdf file
    print('writing file')
    for varname in varkey_list:
        data_xr_var = data_xr[varname]
        data_xr_var.attrs['standard_name'] = var_dict[varname]['standard_name'] #TODO: original long_name and units attrs are now conserved, so do not need to be enforced
        data_xr_var.attrs['coordinates'] = 'longitude latitude'
        #data_xr_var.encoding['zlib'] = True #TODO: way smaller, but also way slower file writing, also with file reading? >> rescale+int16 is better?
        filename = f'ERA5_CDS_atm_{varname}_{dt.datetime.strftime(date_start_zero, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc'
        file_out = os.path.join(dir_output, filename)
        data_xr_var.to_netcdf(file_out)


if __name__ == "__main__":
    if len(os.sys.argv)>1:
        yr=os.sys.argv[1]
        input_dir=os.sys.argv[2]        
    else:
        # yr = '1960'
        # input_dir = './TEMP_meteo_ERA5'
        raise RuntimeError('No arguments were provided\nFirst argument should indicate year as "yyyy".\n Second argument for input_dir') #TODO: re-enable this error again
    convert2FM(yr,input_dir)
