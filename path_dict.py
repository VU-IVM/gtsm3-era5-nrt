# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 13:34:37 2023

@author: veenstra
"""


path_dict_snellius = {'meteo_raw':    '/gpfs/work1/0/einf3499/meteo_ERA5_extended',
                      'meteo_fm':     '/gpfs/work1/0/einf3499/meteo_ERA5_extended_fm',
                      'tides_CDS':    '/gpfs/work1/0/einf3499/tides_CDS_extended',
                      'modeltemplate':'/gpfs/work1/0/einf3499/model_runs_extended/model_template/', # should contain model_input_template and model_files folders #TODO: was model_runs_extended/model_template/model_input_template/ and model_runs/model_template/model_files/'
                      'modelruns':    '/gpfs/work1/0/einf3499/model_runs/slr_tide_surge_runs',
                      'meteo_SLR':    '/gpfs/work1/0/einf3499/meteo_SLR',
                      'meteo_msl':    '/gpfs/work1/0/einf3499/meteo_msl',
                      }

path_dict_h6 =       {'meteo_raw':    './meteo_ERA5',
                      'meteo_fm':     './meteo_ERA5_fm',
                      'tides_CDS':    './tides_CDS',
                      'modeltemplate':'./model_runs_extended/model_template/', # should contain model_input_template and model_files folders #TODO: was model_runs_extended/model_template/model_input_template/ and model_runs/model_template/model_files/'
                      'modelruns':    './model_runs/slr_tide_surge_runs',
                      'meteo_SLR':    './meteo_SLR', #TODO: missing
                      'meteo_msl':    './meteo_msl', #TODO: missing
                      }

path_dict = path_dict_h6
