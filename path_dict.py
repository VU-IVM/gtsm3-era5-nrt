# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:34:37 2023

@author: veenstra
"""

base_dir = '/gpfs/work1/0/einf3499' #snellius

path_dict = {'git_checkout': base_dir+'/00_scripts_git', 
             'meteo_raw':    base_dir+'/01_meteo_ERA5', 
             'meteo_fm':     base_dir+'/02_meteo_ERA5_FM', 
             'tides_CDS':    base_dir+'/03_tides_CDS', 
             'modeltemplate':base_dir+'/00_scripts_git/model_input_template',
             'modelfiles':   base_dir+'/06_model_runs/model_input/model_input_files',
             'modelruns':    base_dir+'/06_model_runs/slr_tide_surge_runs',
             'meteo_SLR':    base_dir+'/05_meteo_SLR',
             'meteo_msl':    base_dir+'/04_meteo_msl',
             }

