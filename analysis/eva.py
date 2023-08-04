import numpy as np
import pandas as pd
from scipy.stats import genpareto as gepd
idx = pd.IndexSlice
import xarray as xr
import os
from gesla import GeslaDataset
from pathlib import Path

def read_gesla(idir: Path):
    '''source: https://gesla787883612.wordpress.com/downloads/ [Aug 2023]  ''' 
    meta_file = os.path.join(idir,"GESLA3_ALL.csv")
    data_path = os.path.join(idir,"GESLA3.0_ALL//")
    filenames = os.listdir(data_path)
    obj_gesla = GeslaDataset(meta_file=meta_file, data_path=data_path)
    return obj_gesla, filenames
    
def process_gesla(obj_gesla, ifile: Path, settings: dict(), message=False):
    data, metadata = obj_gesla.file_to_pandas(ifile)
    #check metadata
    if message==True:
        print('-------------------------------------------------')
        print('Processing station: ' + metadata.site_name)
        print('Overall record quality: '+metadata.overall_record_quality)
        print('Total years of observations: ', metadata.number_of_years ) 
        print('Gauge type: ' + metadata.gauge_type)      
    select = any([metadata.overall_record_quality!= 'No obvious issues', metadata.number_of_years < settings['record_length'],metadata.gauge_type != 'Coastal'])
    gdata = xr.Dataset()
    if select == False:
        # filter data and resample to hourly
        data = data[data['use_flag']==1] 
        data = data.loc['{}-1-1'.format(settings['yearmin']):'{}-1-1'.format(settings['yearmax'])]
        data['time']=data.index 
        data_hourly = data.resample('H',on='time').max()
        data_hourly = data_hourly.dropna(subset='sea_level',how='all')
        data_hourly['dt'] = pd.to_datetime(data_hourly.index).strftime("%Y")
        # count missing data per year
        missing_data = data_hourly.sea_level.groupby(data_hourly.dt).transform('count')
        data_hourly = data_hourly[(missing_data/(365.25*24))*100 > settings['prc_data_missing']]         
        length = len(data_hourly.sea_level) < settings['record_length']*365.25*24
        if length==False: 
            #  convert to xarray 
            data = data_hourly.to_xarray()
            data = data.assign_coords({"station_x_coordinate": metadata.longitude})
            data = data.assign_coords({"station_y_coordinate": metadata.latitude})
            data = data.assign_coords({"station_name": metadata.site_name})
            gdata = data.sea_level.to_dataset()   
    return gdata

def stats(ds_gesla: xr.Dataset ,prcts):
  ds_gesla_stats = ds_gesla.sea_level_detrended.quantile(prcts, dim=('time')) 
  #ds_mean = ds_gesla.sea_level_detrended.mean(dim=('time'))
  #ds_std = ds_gesla.sea_level_detrended.std(dim=('time'))
  #ds_vmin = ds_gesla.sea_level_detrended.min(dim=('time'))
  #ds_vmax = ds_gesla.sea_level_detrended.min(dim=('time'))
  #ds_count = ds_gesla.sea_level_detrended.count(dim=('time'))
  return ds_gesla_stats
    
def detrend(ds_gesla: xr.DataArray, plot = False):
  ''' remove annual means and overall mean '''
  ds_gesla = ds_gesla.assign_coords(year=ds_gesla.time.dt.strftime("%Y"))
  ds = (ds_gesla.groupby("year") - ds_gesla.groupby("year").mean("time"))
  ds_gesla['sea_level_detrended'] = ds['sea_level'] - ds['sea_level'].mean()
  if plot == True:
      fig, axs = plt.subplots(nrows=2)
      ds_gesla.sea_level.plot.line(x='time',ax=axs[0], add_legend=False)   
      ds_gesla.sea_level_detrended.plot.line(x='time',ax=axs[1],add_legend=False)  
  return ds_gesla

def compute_eva(var,istation):
    var = var.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
    probY, potmodel = peak_over_threshold(var)
    potmodel = potmodel.append(pd.DataFrame([istation],index=['station'],columns=potmodel.columns))
    potmodel = potmodel.transpose().set_index('station')
    return potmodel

def pot(var,indep=24*3,res=1,resunit='hour'):
    #    dt = pd.infer_freq(var.index) # way to get the frequency of the data
    #    indep = 4*24# independency test in hours if data is hourly
    #    mindur = 6h*60min minimun duration in units of var above the threshold for the event to be considered a storm
    #    res = 10 min resolution of the data in minutes
    #    resunit = 'min' resolution of the data units
    res = 1
    resunit = 'hour' 
    q = 0.99 # quantile for estimating the threshold
    indep = 3 * 24
        
    #1 - Obtain storms (going up and down the threshold)   
    varth = var.quantile(q) # threshold
    vare = var[var>=varth]
    istorm = np.insert(np.where((vare.index[1:]-vare.index[:-1])>pd.Timedelta(res,resunit))[0],0,-1)

    indx=[];vmax=[];dur=[];
    for i,j in zip(istorm[:-1]+1,istorm[1:]):
        indx.append(vare[i:j+1].idxmax())
        vmax.append(vare[i:j+1].max())
        dur.append(vare[i:j+1].count())
    indx = pd.Index(indx)
    if len(indx)==0:
        print('no INDX', flush=True)
        return np.nan,np.nan,np.nan # stops the computations if the data is wrong at this node
    #2 Check for independence
    indx2=[indx[0]];vmax2=[vmax[0]];dur2=[dur[0]]
    for i,v,d in zip(indx[1:],vmax[1:],dur[1:]):
        if (i-indx2[-1])>pd.Timedelta(indep,'h'):
            indx2.append(i);vmax2.append(v);dur2.append(d)
        else:
            imx = np.argmax([vmax2[-1],v])
            indx2[-1] = [indx2[-1],i][imx]
            vmax2[-1] = [vmax2[-1],v][imx]
            dur2[-1] = dur2[-1]+d
    vare = pd.DataFrame(index=indx2,data=np.vstack((vmax2,dur2)).T,columns=['max','duration'])
    events = vare['max'].groupby(vare.index.year).agg('count')
    return vare,events,varth   

def peak_over_threshold(var):
    ### block_maxima ###
    # Analysis of extremes following block maxima
    # ds: pandas dataframe
    # rv: return values to compute the modelled return periods
    # fitd: distribution function to use for the extreme value analysis
    # qopt: type of threshold indicated, value or quantile based
    # indep: independence between events [see definition of pot for the units]
    # nr: set of bootstrap samples of size ns
    # per: percentiles to return parameter estimations
    rv = np.array([1,2,5,10,25,50,100])
    q = 0.99
    # Empirical model
    #1- find storms
    potmax,evts,thresh = pot(var) 
    if np.isnan(thresh):
        return pd.DataFrame(),pd.DataFrame() # stops the computations if the data is wrong at this node
    
    #2- sort data smallest to largest    
    pmord = potmax.sort_values(by='max');pmord.index.name = 'date'
    pmord = pmord.reset_index()
    pmord = pmord.rename(columns={'max':'emp'})

    #3- count total obervations
    n = len(pmord.index)
    pmord['probex'] = ((n - (pmord.index.values+1) + 1) / (n + 1))          
    pmord['rv'] = 1./evts.mean()*(1./pmord['probex'])   #1/lambda*1/F lambda is Poisson parameter the average number of events per year  

    #4- Fit gpd   
    dpars = ['lambda','shape','location','scale', 'thresh', 'qthres']
    dparsrv = np.hstack((dpars,rv.astype(str)))
    pgpd = np.array(dparsrv)
       
    potmodel = pd.DataFrame(columns=['mean'],index=pgpd)
    potmodelc = pd.DataFrame(columns=['mean'],index=np.hstack(([rv[0]],np.arange(rv[1],rv.max()+1))))
    potmodel.loc[['shape','location','scale'],'mean'] = gepd.fit(pmord['emp'].values,0,loc=var.quantile(q))
    potmodel.loc[['lambda'],'mean'] = evts.mean()
    potmodel.loc[['thresh'],'mean'] = thresh
    potmodel.loc[['qthres'],'mean'] = var.quantile(q)
    
    potmodelc.loc[:,'mean'] = gepd.ppf(1.-1./(potmodel.loc['lambda','mean']*potmodelc.index.values),potmodel.loc['shape','mean'],loc=potmodel.loc['location','mean'],scale=potmodel.loc['scale','mean'])
    potmodel.loc[rv.astype(str),'mean']=potmodelc.loc[rv,'mean'].values
    pmord['mod_mean'] = gepd.ppf(1.-pmord['probex'].values,potmodel.loc['shape','mean'],loc=potmodel.loc['location','mean'],scale=potmodel.loc['scale','mean'])
    return pmord.drop(columns='probex'),potmodel