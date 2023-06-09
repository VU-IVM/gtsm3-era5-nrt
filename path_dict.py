# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:34:37 2023

@author: veenstra
"""

import sys

# base_dir = '/gpfs/work1/0/einf3499' #snellius
base_dir = '/p/1230882-emodnet_hrsm/GTSMv3.0EMODnet/sealevelmonitor/workflow_test' #h6
if sys.platform=='win32': #on WCF
    base_dir = base_dir.replace('/p','p:')

path_dict = {'meteo_raw':    base_dir+'/meteo_ERA5', #TODO: contains _extended on snellius
             'meteo_fm':     base_dir+'/meteo_ERA5_fm', #TODO: contains _extended on snellius
             'tides_CDS':    base_dir+'/tides_CDS', #TODO: contains _extended on snellius
             'modeltemplate':base_dir+'/model_runs/model_template', # should contain model_input_template and model_files folders #TODO: was model_runs_extended/model_template/model_input_template/ and model_runs/model_template/model_files/' on snellius
             'modelruns':    base_dir+'/model_runs/slr_tide_surge_runs',
             'meteo_SLR':    base_dir+'/meteo_SLR', #TODO: missing on h6
             'meteo_msl':    base_dir+'/meteo_msl', #TODO: missing on h6
             }

