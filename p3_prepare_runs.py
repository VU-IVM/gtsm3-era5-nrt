#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 17:25:23 2018

@author: muis
"""
import os
import sys
import click
from datetime import date, datetime
import numpy as np
import os
import shutil
import fnmatch
import datetime
import templates
import netCDF4 
from print_fou import print_fou
from distutils.dir_util import copy_tree

sys.path.insert(0, '/gpfs/work1/0/einf3499/model_runs/model_template/')
#folder for current script
base_path=os.path.dirname(os.path.realpath(os.sys.argv[0]))
  
@click.command()
@click.option('--output_dir', required=True, type=str,
              help='Output directory for downloaded files',)
@click.option('--date_string', required=True, type=str,
              help='String with year and month of current month (in YYYY_MM format)',)

def run_GTSM_monthly_runs(base_path,date_string):  
  tdate = datetime.strptime(date_string, '%Y_%m').date()
  base_name="ERA5"
  var ='slr_tide_surge'
  
  meteodir="/projects/0/einf3499/meteo_ERA5_fm/" 
  SLRfile="/projects/0/einf3499/meteo_SLR/TotalSeaLevel_MapsSROCC_rcp85_Perc50_zero1986to2005_dflow_extrap.nc"
  MSLdir="/projects/0/einf3499/meteo_msl/" 
  
  templatedir="/gpfs/work1/0/einf3499/model_runs/model_template/model_input_template"
  modelfilesdir="/gpfs/work1/0/einf3499/model_runs/model_template/model_files"
      
  meteofile_u=os.path.join(meteodir,'ERA5_CDS_atm_u10_%s-01-01_%s-01-01.nc' %(int(yr), int(yr)+1))
  meteofile_v=os.path.join(meteodir,'ERA5_CDS_atm_v10_%s-01-01_%s-01-01.nc' %(int(yr), int(yr)+1))
  meteofile_p=os.path.join(meteodir,'ERA5_CDS_atm_msl_%s-01-01_%s-01-01.nc' %(int(yr), int(yr)+1))
        
  MSLcorrfile=os.path.join(MSLdir, 'ERAInterim_average_msl_neg_%s1215_%s0101.nc' %(int(yr-1), int(yr)+1))
      
  run_dir="/gpfs/work1/0/einf3499/model_runs/%s_runs/model_input_%s_%d" %(var,base_name,yr)
      
  # copy model files and template
  #-------------------
  try:
     os.stat(newdir)               
     raise Exception("Directory already exists"+newdir) #this is to avoid overwritting runs that are there already by mistake
  except OSError:
      print("copying ",templatedir," to ",newdir)
      shutil.copytree(templatedir,newdir,symlinks=False,ignore=None)

  copy_tree(modelfilesdir, newdir)

  # calculate start and end times based on chosen reference time
  #-------------------
  refdate=datetime(1900,1,1) 
  spinup=15 #in days
  inidate=(y-timedelta(days=spinup))    
  enddate=years[iyear+1]
  tstart=(inidate-refdate).total_seconds() 
  tend=(enddate-refdate).total_seconds()

  print("reference date is ",str(trefsim))
  print("tstart is ",str(tstart)," hours since ref time") 
  print("tstop  is ",str(tstop)," hours since ref time")         
      
  ## change templates
  #-------------------
  keywords_MDU={'REFDATE':user_refdate.strftime("%Y%m%d"),'TSTART':str(tstart),'TSTOP':str(tend)}
  templates.replace_all(os.path.join(newdir,"gtsm_fine.mdu.template"), os.path.join(newdir,"gtsm_fine.mdu"),keywords_MDU,'%')

  keywords_EXT={'METEOFILE_GLOB_U':meteofile_u,'METEOFILE_GLOB_V':meteofile_v,'METEOFILE_GLOB_P':meteofile_p, 'METEOFILE_MSLCORR':MSLcorrfile,'METEOFILE_SLR':SLRfile}
  templates.replace_all(os.path.join(templates_path,'gtsm_forcing.ext.template'),os.path.join(rundir,extfile),keywords_EXT,'%') 

  keywords_QSUB={'NDOM':str(ndom),'NNODE':str(nnode),'JOBNAME':workfolder} 
  templates.replace_all(os.path.join(templates_path,'%s.template'%(shfile)),os.path.join(rundir,'%s'%(shfile)),keywords_QSUB,'%')

  os.system("cd "+newdir+"; chmod -R 777 *")    
    
if __name__ == "__main__":
  prepare_GTSM_monthly_runs()
