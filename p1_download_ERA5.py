#!/usr/bin/env python

def downdload_era5(startdate,enddate,outdir):
    # import modules
    import os
    import os.path
    import cdsapi
    import datetime
    import netCDF4 as nc
    # download the data
    c=cdsapi.Client()
    startyear=startdate.year
    endyear=enddate.year
    for year in range(startyear,maxyear+1):
        print ('YEAR ',year)
        startmonth=startdate.year
        endmonth=enddate.year
        for month in range(startmonth,endmonth+1):
            lastday=calendar.monthrange(year,month)[1]
            bdate="%s%02d01"%(year,month)
            edate="%s%02d%s"%(year,month,lastday)
            targetfile=os.path.join(outdir,"ERA5_CDS_atm_%s_%s.nc" %(bdate,edate))
            dayarray=["{:02}".format(iday) for iday in range(1,lastday+1)]           
            print ("######### ERA-5 from CDS  #########")
            print ('get data from ', bdate,' to ',edate,' (YYYYMMDD)')
            print ("################################")
            
            # Atm variables
            c.retrieve(
	        'reanalysis-era5-single-levels',
	        {
	            'product_type':'reanalysis',
	            'format':'netcdf',
	            'variable':[
	                '10m_u_component_of_wind','10m_v_component_of_wind','mean_sea_level_pressure'
	            ],
	            'year':str(year),
	            'month':["{:02}".format(month)],
	            'day':dayarray,
	            'time':[
	                '00:00','01:00','02:00',
	                '03:00','04:00','05:00',
	                '06:00','07:00','08:00',
	                '09:00','10:00','11:00',
	                '12:00','13:00','14:00',
	                '15:00','16:00','17:00',
	                '18:00','19:00','20:00',
	                '21:00','22:00','23:00'
	            ]
	        },
           targetfile)
           # change attributes
           ds = xr.open_dataset(targetfile)
           ds.uas.attrs['standard_name'] = 'eastward_wind'
           ds.vas.attrs['standard_name'] = 'northward_wind'


if __name__ == "__main__":
    # read input arguments
    if len(os.sys.argv)>1:
        startdate=datetime.datetime.strptime(os.sys.argv[1], "%Y-%m-%d")
        enddate=datetime.datetime.strptime(os.sys.argv[2], "%Y-%m-%d")
        outdir=os.sys.argv[3])        
    else:
        raise RuntimeError('No arguments were provided\nFirst argument should indicate startdate as "%Y-%m-%d".\n Second argument for enddate. \n Third argument the outdir')
    downdload_era5(startdate,enddate,outdir)