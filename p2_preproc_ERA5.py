#!/usr/bin/env python

import os
import xarray as xr
import datetime as dt
import numpy as np
import glob
from path_dict import path_dict


def recompute_scalefactor_offset_filelist(file_nc_list, data_xr):
    """
    First computes min/max from scalefactor/offset per variable in source ncfiles and keeps track of the extremes.
    Then computes xcalefactor/offset per variable in merged ncfile and overwrites the encoding.
    This makes sure int16 variables can be written while making sure the entire range is represented.
    """
    def compute_scale_and_offset(min, max, n=16):
        # stretch/compress data to the available packed range
        scale_factor = (max - min) / (2 ** n - 1)
        # translate the range to be symmetric about zero
        add_offset = min + 2 ** (n - 1) * scale_factor
        return (scale_factor, add_offset)
    #scale_factor, add_offset = compute_scale_and_offset(-100, 100)
    
    def compute_min_and_max(scale_factor,add_offset,n=16):
        vals_range = scale_factor * (2 ** n - 1)
        min = add_offset - 2 ** (n - 1) * scale_factor
        max = min + vals_range
        return min, max
    #vals_min, vals_max = compute_min_and_max(scale_factor,add_offset)
    
    varkey_list = list(data_xr.data_vars)
    
    #get lowest min and highest max from all files for all variables
    min_dict = {x:np.nan for x in varkey_list}
    max_dict = {x:np.nan for x in varkey_list}
    for file_nc in file_nc_list:
        ds = xr.open_dataset(file_nc)
        for varn in varkey_list:
            nbits = ds[varn].encoding['dtype'].itemsize*8 #itemsize is in bytes, multiply by 8 to get bits
            var_min,var_max = compute_min_and_max(ds[varn].encoding['scale_factor'], ds[varn].encoding['add_offset'], n=nbits)
            min_dict[varn] = np.nanmin([min_dict[varn],var_min])
            max_dict[varn] = np.nanmax([max_dict[varn],var_max])
    
    #recompute scale_factor and add_offset for all variables
    for varn in varkey_list:
        nbits = data_xr[varn].encoding['dtype'].itemsize*8 #itemsize is in bytes, multiply by 8 to get bits
        var_sf, var_offset = compute_scale_and_offset(min_dict[varn],max_dict[varn],n=nbits)
        #print(varn, var_sf, var_offset)
        print(f"replacing {varn} scale_factor: {data_xr[varn].encoding['scale_factor']} with {var_sf}")
        print(f"replacing {varn} add_offset: {data_xr[varn].encoding['add_offset']} with {var_offset}")
        
        data_xr[varn].encoding['scale_factor'] = var_sf
        data_xr[varn].encoding['add_offset'] = var_offset


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
    dir_output = input_dir.replace('meteo_ERA5','meteo_ERA5_fm')
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    
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
    
    recompute_scalefactor_offset_filelist(file_nc_list, data_xr)
    
    
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
    # if time_slice.stop not in times_pd.index: #TODO: re-enable this part of the code
    #     raise Exception(f'ERROR: time_slice_stop="{time_slice.stop}" not in selected files, timerange: "{times_pd.index[0]}" to "{times_pd.index[-1]}"')
    
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
    #TODO: doing this drops all encoding from variables, causing them to be converted into floats. Also makes sense since 0 pressure does not fit into int16 range as defined by scalefac and offset
    #'scale_factor': 0.17408786412952254, 'add_offset': 99637.53795606793
    #99637.53795606793 - 0.17408786412952254*32768
    #99637.53795606793 + 0.17408786412952254*32767
    if zerostart: #TODO: is one zero field not enough? (is replacing first field not also ok? (results in 1hr transition period)
        field_zerostart = xr.zeros_like(data_xr.isel(time=[0,0])) #two times first field with zeros, encoding of data_vars is dropped)
        field_zerostart['time'] = [times_pd.index[0]-dt.timedelta(days=2),times_pd.index[0]-dt.timedelta(days=1)] #this drops time var encoding
        #for varn in varkey_list+['time']: #loop over variables and overwrite encoding with original #TODO: 0-msl values do not fit in int32 range, so will be incorrect
        #    field_zerostart[varn].encoding = data_xr[varn].encoding
        data_xr = xr.concat([field_zerostart,data_xr],dim='time',combine_attrs='no_conflicts') #combine_attrs argument prevents attrs from being dropped
    
    # create variables 'dictionary' with attrs
    var_dict = {
      "u10" : {
        "standard_name" : "eastward_wind", #TODO: missing from ERA5 dataset >> couple via quantity/varname in extfile
        "scale_factor" : float(0.01), #TODO: scale_factor and add_offset are not written to file, should also be written towards encoding instead of attrs. Also dtype and _FillValue encoding should be (over)written
        "add_offset" : float(0)},
      "v10" : {
        "standard_name" : "northward_wind",#TODO: missing from ERA5 dataset
        "scale_factor" : float(0.01),
        "add_offset" : float(0)},
      "msl" : {
        "standard_name" : "air_pressure", #TODO: is air_pressure_at_mean_sea_level in ERA5 dataset
        "scale_factor" : float(1),
        "add_offset" : float(100000)}} #TODO: zero-fields will not fit in the resulting int16 range, so dtype=int16 is not possible
    #write to netcdf file
    print('writing file')
    for varname in varkey_list[:1]: #TODO:revert back to entire list
        data_xr_var = data_xr[varname]
        data_xr_var.attrs['standard_name'] = var_dict[varname]['standard_name'] #TODO: original long_name and units attrs are now conserved, so do not need to be enforced
        data_xr_var.attrs['coordinates'] = 'longitude latitude'
        #TODO: reference speed/disksize for one float32 variable: 2min, 4.6GB
        #data_xr_var.encoding['zlib'] = True #TODO: way smaller, but also way slower file writing (with complevel>0), also with file reading?
        #data_xr_var.encoding['compression'] = "gzip"
        # if varname in ['u10','v10']: # speed/disksize for one int16 variable: 1min, 2.3GB
        #     data_xr_var.encoding['scale_factor'] = 0.01 # TODO: scale_factor=0.01 this results in an accuracy of 1cm/s in x an y direction. Rounding the separate components might slightly shift the wind direction
        #     data_xr_var.encoding['add_offset'] = 0
        #     data_xr_var.encoding['dtype'] = 'int16'
        #     data_xr_var.encoding['_FillValue'] = -32767
        filename = f'ERA5_CDS_atm_{varname}_{dt.datetime.strftime(date_start_zero, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc'
        file_out = os.path.join(dir_output, filename)
        data_xr_var.to_netcdf(file_out)
        
        ds_out = xr.open_dataset(file_out)
        import matplotlib.pyplot as plt
        fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(12,8))
        for iA,ax in enumerate([ax1,ax2,ax3]):
            ds_out.msl.isel(time=iA).plot(ax=ax)


if __name__ == "__main__":
    if len(os.sys.argv)>1:
        yr=os.sys.argv[1]
    else:
        yr = '1960'
        #raise RuntimeError('No arguments were provided\nFirst argument should indicate year as "yyyy".')
    tstart = dt.datetime.now()
    convert2FM(yr)
    print('time passed:',dt.datetime.now()-tstart)
    
