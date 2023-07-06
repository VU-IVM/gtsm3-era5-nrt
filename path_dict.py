# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:34:37 2023

@author: veenstra
"""

import sys

base_dir = '/gpfs/work1/0/einf3499' #snellius
#base_dir = '/p/1230882-emodnet_hrsm/GTSMv3.0EMODnet/sealevelmonitor/workflow_test' #h6
if sys.platform=='win32': #on WCF
    base_dir = base_dir.replace('/p','p:')

path_dict = {'base_dir':     base_dir,
             'meteo_raw':    base_dir+'/01_meteo_ERA5', 
             'meteo_fm':     base_dir+'/02_meteo_ERA5_FM', 
             'tides_CDS':    base_dir+'/03_tides_CDS', 
             'modeltemplate':base_dir+'/06_model_runs/model_input', # should contain model_input_template and model_files folders 
             'modelruns':    base_dir+'/06_model_runs/slr_tide_surge_runs',
             'meteo_SLR':    base_dir+'/05_meteo_SLR', #TODO: missing on h6
             'meteo_msl':    base_dir+'/04_meteo_msl', #TODO: missing on h6
             }

