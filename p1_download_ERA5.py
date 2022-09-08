#!/usr/bin/env python

def download_era5(tstart,outdir):
  #tstart='01-01-2021'
  # import modules
  import os
  import os.path
  from datetime import date, timedelta, datetime
  import calendar
  import cdsapi
  import netCDF4 as nc
  # find times for monthly downlods
  tstart = datetime.strptime(tstart, "%d-%m-%Y").date()
  lastday=calendar.monthrange(tstart.year,tstart.month)[1]
  dayarray=["{:02}".format(iday) for iday in range(1,lastday+1)]     
  # daily download
  print ("######### ERA-5 from CDS  #########")
  print ('get data from ', tstart,'(1 Day)')
  print ("################################")
  # I/O - download the data
  for day in dayarray:
    print(day)
    tday = datetime(day=int(day),month=tstart.month,year=tstart.year).date()
    targetfile=os.path.join(outdir,"ERA5_CDS_atm_%s.nc" %(tday))
    if os.path.isfile(targetfile)==False: 
      c=cdsapi.Client()
      c.retrieve('reanalysis-era5-single-levels',
          {'product_type':'reanalysis',
          'format':'netcdf',
          'variable':['10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure'],
          'year':tstart.year,
          'month':["{:02}".format(tstart.month)],
          'day':day,
          'time':['00:00','01:00','02:00','03:00','04:00','05:00',
                  '06:00','07:00','08:00','09:00','10:00','11:00',
        	        '12:00','13:00','14:00','15:00','16:00','17:00',
                  '18:00','19:00','20:00','21:00','22:00','23:00']
          },targetfile)
   
def convert2FM(ifile,ofile=None):
    """
    This function converts a ERA5 file from CDS to FM format
    """
    # import modules
    import xarray as xr
    import datetime as dt
    import numpy as np
    import pandas as pd
    # import dataset
    ds = xr.open_dataset(ifile)

    # copy latitude and create new longitude
    lats = ds['latitude'][:].values
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
    # v10
    datasets = []
    for i in range(0,len(rel_times)):
      v10 = ds.v10.isel(time=i).values
      var_new=np.hstack((v10[:,part1],v10[:,part2]))
      coords =  {'lat': lats, 'lon': lons_new, 'time': rel_times[0]}
      da = xr.DataArray(var_new, coords=coords, dims=['lat', 'lon'])
      da.name = 'v10'
      dat = xr.concat([da], 'time')
      datasets.append(dat)
    ds_v10_merged = xr.concat(datasets, dim='time') 
    # u10
    datasets = []
    for i in range(0,len(rel_times)):
      u10 = ds.u10.isel(time=i).values
      var_new=np.hstack((u10[:,part1],u10[:,part2]))
      coords =  {'lat': lats, 'lon': lons_new, 'time': rel_times[0]}
      da = xr.DataArray(var_new, coords=coords, dims=['lat', 'lon'])
      da.name = 'u10'
      dat = xr.concat([da], 'time')
      datasets.append(dat)
    ds_u10_merged = xr.concat(datasets, dim='time') 
    # mslp
    datasets = []
    for i in range(0,len(rel_times)):
      msl = ds.msl.isel(time=i).values
      var_new=np.hstack((msl[:,part1],msl[:,part2]))
      coords =  {'lat': lats, 'lon': lons_new, 'time': rel_times[0]}
      da = xr.DataArray(var_new, coords=coords, dims=['lat', 'lon'])
      da.name = 'msl'
      dat = xr.concat([da], 'time')
      datasets.append(dat)
    ds_msl_merged = xr.concat(datasets, dim='time') 
    # merge and set attributes
    ds_all = xr.merge([ds_v10_merged,ds_u10_merged,ds_msl_merged]) 
    print(ds_all)
    print(ds_all.u10)
    print(ds.u10)
# to do
# attributes, scale + offset
        # change attributes    
    ds_all.u10.attrs['standard_name'] = 'eastward_wind'
    ds_all.v10.attrs['standard_name'] = 'northward_wind'
    ds_all.msl.attrs['standard_name'] = 'air_pressure'
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
    if ofile==None:
      ofile = ifile.replace('.nc','_fm.nc')
    ds_all.to_netcdf(ofile.replace, encoding=encoding)
#        overlapleft=(x<-179) #overlap
#        x_new=np.hstack((x[overlapright]-360,x,x[overlapleft]+360))
#        nx=len(x_new)
#        longitude_dim=dest.createDim('longitude',nx)
#        longitude_type=maps.getVariableType('longitude')
#        new_longitude=dest.createVariable('longitude',longitude_type,('longitude'))
#        dest.copyVariableAttributesFrom(maps,'longitude')
#        new_longitude[:]=x_new;
#    
#    
#    if os.path.isfile(ofile):
#      raise RuntimeError('Output file already exists. Please remove and try again. File name is '+newfilename)
#      #os.remove(newfilename)
#
#    maps=EraMaps(filename,'r')
#    times=maps.getRelativeTimes()
#    print "Postprocessing " + filename
#    first=True
#    tcount=0
#    if first:
#        dest=EraMaps(newfilename,'NETCDF4','w')
#        dest.copyGlobalDataFrom(maps,skipDims=['longitude','time'])
#        dest.copyVariableFrom(maps,'latitude')
#        #create Dimesions and variables
#        #longitude
#
#        
#        #time
#        time_dim = dest.createDim('time',None)
#        time_type=maps.getVariableType('time')
#        new_time=dest.createVariable('time',time_type,('time'))
#        dest.copyVariableAttributesFrom(maps,'time')
#        dest.createVariableAttribute('time','standard_name','time')
#        
#        #msl
#        msl_type=maps.getVariableType('msl')
#        dest.createVariable('msl',msl_type,('time','latitude','longitude'))
#        dest.copyVariableAttributesFrom(maps,'msl')
#        dest.createVariableAttribute('msl','standard_name','air_pressure')
#        dest.createVariableAttribute('msl','coordinates','lon lat')
#        dest.createVariableAttribute('msl','scale_factor',scalef_msl)
#        dest.createVariableAttribute('msl','add_offset',offset_msl)
#        
#        #u10
#        u10_type=maps.getVariableType('u10')
#        dest.createVariable('u10',u10_type,('time','latitude','longitude'))
#        dest.copyVariableAttributesFrom(maps,'u10')
#        dest.createVariableAttribute('u10','standard_name','eastward_wind')
#        dest.createVariableAttribute('u10','coordinates','lon lat')
#        dest.createVariableAttribute('u10','scale_factor',scalef_u10)
#        dest.createVariableAttribute('u10','add_offset',offset_u10)
#        #v10
#        v10_type=maps.getVariableType('v10')
#        dest.createVariable('v10',u10_type,('time','latitude','longitude'))
#        dest.copyVariableAttributesFrom(maps,'v10')
#        dest.createVariableAttribute('v10','standard_name','northward_wind')
#        dest.createVariableAttribute('v10','coordinates','lon lat')
#        dest.createVariableAttribute('v10','scale_factor',scalef_v10)
#        dest.createVariableAttribute('v10','add_offset',offset_v10)
#        first=False
#        
#    # copy data
#    ntimes=len(times)
#    dest.setRelativeTimes(times,tcount)    
#  
#    for itime in np.arange(ntimes):
#        pressure=maps.getPressure(itime)
#        #move 180:360 part to -180:0 so field now runs from longitute -180 to 180
#        new_pressure=np.hstack((pressure[:,overlapright],pressure[:,:],pressure[:,overlapleft]))
#        dest.setPressure(tcount,new_pressure)
#        
#        xwind=maps.getXVelocity(itime)
#        new_xwind=np.hstack((xwind[:,overlapright],xwind[:,:],xwind[:,overlapleft]))
#        dest.setXVelocity(tcount,new_xwind)
#        
#        ywind=maps.getYVelocity(itime)
#        new_ywind=np.hstack((ywind[:,overlapright],ywind[:,:],ywind[:,overlapleft]))
#        dest.setYVelocity(tcount,new_ywind)
#        tcount=tcount + 1
#            
#        #clean up
#    maps.close()
#dest.close()
#
#
#ds.to_netcdf(outfile, encoding={"my_var": {
#    "dtype": 'int16',
#    "scale_factor": scale_factor,
#    "add_offset": add_offset,
#    "_FillValue": -32767,
#}})

if __name__ == "__main__":
    # read input arguments
    import os
    import datetime
    if len(os.sys.argv)>0:
      tstart=os.sys.argv[1]
      outdir=os.sys.argv[2]        
    else:
      raise RuntimeError('No arguments were provided\nFirst argument should indicate startdate as "%Y-%m-%d".\n Second argument for outdir. Script will download monthly files per day')
    download_era5(tstart,outdir)
    convert2FM('ERA5_CDS_atm_2015-01-01.nc')