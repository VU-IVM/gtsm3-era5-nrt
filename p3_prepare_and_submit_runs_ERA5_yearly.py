#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 17:25:23 2018

@author: muis
"""
def run_GTSM_yearly_runs(yr0,yr1):
  # import modules
  import sys
  sys.path.insert(0, '/gpfs/work1/0/einf3499/model_runs/model_template/')
  import numpy as np
  import os
  import shutil
  import fnmatch
  import datetime
  import templates
  import netCDF4 
  from print_fou import print_fou
  from distutils.dir_util import copy_tree
  import glob 
  #folder for current script
  base_path=os.path.dirname(os.path.realpath(os.sys.argv[0]))
  
  #USER SETTINGS
  yrs = np.arange(yr0,yr1+1)
  base_name="ERA5"
  
  #PATHS
  meteodir="/projects/0/einf3499/meteo_ERA5_fm/" 
  SLRfile="/projects/0/einf3499/meteo_SLR/TotalSeaLevel_MapsSROCC_rcp85_Perc50_zero1986to2005_dflow_extrap.nc"
  MSLdir="/projects/0/einf3499/meteo_msl/" 
  
  templatedir="/gpfs/work1/0/einf3499/model_runs/model_template/model_input_template"
  modelfilesdir="/gpfs/work1/0/einf3499/model_runs/model_template/model_files"
  var ='slr_tide_surge'
  
  trefsim=datetime.datetime(1900,1,1)
  
  # loop over years and submit  
  for iyr in range(len(yrs)):       
      yr = yrs[iyr]     
      meteofile_u=os.path.join(meteodir,'ERA5_CDS_atm_u10_%s-01-01_%s-01-01.nc' %(int(yr), int(yr)+1))
      meteofile_v=os.path.join(meteodir,'ERA5_CDS_atm_v10_%s-01-01_%s-01-01.nc' %(int(yr), int(yr)+1))
      meteofile_p=os.path.join(meteodir,'ERA5_CDS_atm_msl_%s-01-01_%s-01-01.nc' %(int(yr), int(yr)+1))
      print(meteofile_u)
      
      MSLcorrfile=os.path.join(MSLdir, 'ERAInterim_average_msl_neg_%s1215_%s0101.nc' %(int(yr-1), int(yr)+1))
      
      print("   ",base_name)
      newdir="/gpfs/work1/0/einf3499/model_runs/%s_runs/model_input_%s_%d" %(var,base_name,yr)
      
      #create new folder and copy contents of template folder there
      try:
         os.stat(newdir)               
         raise Exception("Directory already exists"+newdir) #this is to avoid overwritting runs that are there already by mistake
      except OSError:
          print("copying ",templatedir," to ",newdir)
          shutil.copytree(templatedir,newdir,symlinks=False,ignore=None)
          
      #copy static model files
      copy_tree(modelfilesdir, newdir)
      
      #calculate start and end times based on chosen reference time
      refdate=datetime(1900,1,1) 
      spinup=15 #in days
      inidate=(y-timedelta(days=spinup))    
      enddate=years[iyear+1]
      tstart=(inidate-refdate).total_seconds() 
      tend=(enddate-refdate).total_seconds()

      print("#######################################")
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
      # start the run
      #os.system("./sbatch_snellius_delft3dfm2022.04_1x128.sh")
    
if __name__ == "__main__":
  # read input arguments
  import os
  if len(os.sys.argv)>0:
    yr0=int(os.sys.argv[1])
    yr1=int(os.sys.argv[2])
    print(yr0,yr1)       
  else:
    raise RuntimeError('No arguments were provided\nFirst argument should indicate startdate as "%Y-%m-%d".\n Second argument for enddate and thirs for output directory. Script will merge daily files and convert to FM format')
  run_GTSM_yearly_runs(yr0,yr1)
