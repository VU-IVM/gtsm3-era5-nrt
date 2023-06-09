#! /usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 17:25:23 2018

@author: muis
"""

import os
import datetime as dt
import shutil
from distutils.dir_util import copy_tree


def prepare_GTSM_monthly_runs(yr,base_dir):
    # calculate start and end times based on chosen reference time
    tdate = dt.datetime.strptime(yr,'%Y').date()
    date_start = dt.datetime(tdate.year,1,1)-dt.timedelta(days=17) # imposed 1 day zero, 1 day transition, 15 days spinup 
    date_end = dt.datetime(tdate.year+1,1,1)
    
    refdate = dt.datetime(1900,1,1)
    tstart = (date_start-refdate).total_seconds()
    tstop = (date_end-refdate).total_seconds()
    
    print("reference date is ",str(refdate))
    print("tstart is ",str(tstart)," seconds since ref time")
    print("tstop  is ",str(tstop)," seconds since ref time")
    # files
    meteofile_u=os.path.join(base_dir,f'meteo_ERA5_fm_extended/ERA5_CDS_atm_u10_{dt.datetime.strftime(date_start, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc')
    meteofile_v=os.path.join(base_dir,f'meteo_ERA5_fm_extended/ERA5_CDS_atm_v10_{dt.datetime.strftime(date_start, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc')
    meteofile_p=os.path.join(base_dir,f'meteo_ERA5_fm_extended/ERA5_CDS_atm_msl_{dt.datetime.strftime(date_start, "%Y-%m-%d")}_{dt.datetime.strftime(date_end, "%Y-%m-%d")}.nc')
    
    SLRfile=os.path.join(base_dir,"meteo_SLR/TotalSeaLevel_MapsSROCC_rcp85_Perc50_zero1986to2005_dflow_extrap.nc")  
    MSLcorrfile=os.path.join(base_dir, "meteo_msl/ERAInterim_average_msl_neg_%s1215_%s0101.nc" %(int(tdate.year-1), int(tdate.year)+1)) # need to be updated
    
    templatedir=os.path.join(base_dir, "model_runs_extended/model_template/model_input_template")
    modelfilesdir=os.path.join(base_dir, "model_runs/model_template/model_files")
    run_dir= os.path.join(base_dir, f"model_runs_extended/slr_tide_surge_runs/model_input_ERA5_{tdate.year}")
    # os.makedirs(run_dir, exist_ok=True)
    
    # copy model files and template
    try:
        os.stat(run_dir)
        raise Exception("Directory already exists ", run_dir)
    except OSError:
        print("copying ",templatedir," to ",run_dir)
        shutil.copytree(templatedir,run_dir,symlinks=False,ignore=None)
    copy_tree(modelfilesdir, run_dir)
    
    # change templates
    keywords_MDU={'%REFDATE%':refdate.strftime("%Y%m%d"),'%TSTART%':str(tstart),'%TSTOP%':str(tstop)}
    replace_all(os.path.join(run_dir,"gtsm_fine.mdu.template"), os.path.join(run_dir,"gtsm_fine.mdu"),keywords_MDU)
    
    keywords_EXT={'%METEOFILE_GLOB_U%':meteofile_u,'%METEOFILE_GLOB_V%':meteofile_v,'%METEOFILE_GLOB_P%':meteofile_p, '%METEOFILE_MSLCORR%':MSLcorrfile,'%METEOFILE_SLR%':SLRfile}
    replace_all(os.path.join(templatedir,'gtsm_fine.ext.template'),os.path.join(run_dir,"gtsm_fine.ext"),keywords_EXT) 
    
    shfile = 'sbatch_snellius_delft3dfm2022.04_1x128cores_yearly.sh'
    workfolder=f"ERA5_{tdate.year}"
    keywords_QSUB={'%JOBNAME%':workfolder}
    replace_all(os.path.join(templatedir,shfile),os.path.join(run_dir,'%s'%(shfile)),keywords_QSUB)
    
    os.system("cd "+run_dir+"; chmod -R 777 *")


def replace_all(template_name,out_name,replace_dict):
    '''replace all occurrances of replace_dict keys with its values'''
    f_in  = open(template_name,'r')
    f_out = open(out_name,'w')
    data = f_in.read() # string of all file content
    for i, j in replace_dict.items(): #replace per key
        data = data.replace(i,j)
    f_out.write(data)
    f_in.close()
    f_out.close()
    
    #check for missed replacements
    msg = []
    if '%' in data:
        for line in data.split('\n'):
            if '%' in line:
                msg.append(line)
        msg_newlines = '\n'.join(msg)
        raise Exception(f'Not all keys replaced:\n{msg_newlines}')


if __name__ == "__main__":
    if len(os.sys.argv)>1:
        yr=os.sys.argv[1]
        base_dir=os.sys.argv[2]        
    else:
        # base_dir = '.'
        # run_dir= os.path.join(base_dir, "template_files")
        # refdate = dt.datetime(2011,1,1)
        # tstart = '0'
        # tstop = '100'
        # keywords_MDU={'%REFDATE%':refdate.strftime("%Y%m%d"),'%TSTART%':str(tstart),'%TSTOP%':str(tstop)}
        # replace_all(os.path.join(run_dir,"gtsm_model.mdu.template"), os.path.join(run_dir,"gtsm_model.mdu"),keywords_MDU)   
        raise RuntimeError('No arguments were provided')
    prepare_GTSM_monthly_runs(yr,base_dir)
