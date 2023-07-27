import numpy as np
import pandas as pd
from scipy.stats import genpareto as gepd
idx = pd.IndexSlice

def compute_eva(var):
    var = var.to_dataframe().loc[:, 'sea_level_detrended'].dropna()
    probY, potmodel = peak_over_threshold(var)
    return probY

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