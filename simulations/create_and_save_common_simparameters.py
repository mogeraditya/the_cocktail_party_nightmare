# -*- coding: utf-8 -*-
"""Create a common parameter set to initialise all simulations and 
change only the parameters of interest for each simulation type.

Created on Sat Jun 08 18:59:27 2019

@author: tbeleyur
"""

import dill 
import numpy as np 
import pandas as pd
import statsmodels.api as sm
import os
print(os.getcwd())
os.chdir("../data/")
bistatic_TS_file = 'bistatic_TS_bat.csv'
shadowing_model_file = 'acoustic_shadowing_model_p3.pkl'
tempmasking_file = 'temporal_masking_function_p3.pkl'
spatial_unmasking_file = 'spatial_release_fn.csv'


kwargs={}
kwargs['interpulse_interval'] = 0.1
kwargs['v_sound'] = 330.0
kwargs['simtime_resolution'] = 10**-6
kwargs['echocall_duration'] = 0.0025        
reflection_func = pd.read_csv(bistatic_TS_file)
kwargs['reflection_function'] = reflection_func
kwargs['heading_variation'] = 10.0
kwargs['min_spacing'] = 0.5
kwargs['Nbats'] = None
kwargs['source_level'] = {'dBSPL' : 100, 'ref_distance':1.0}
kwargs['hearing_threshold'] = 20
kwargs['rectangle_width'] = 0.25
kwargs['implement_shadowing'] = True

# load the model to predict shadowing 

kwargs['acoustic_shadowing_model']=sm.load(shadowing_model_file)


with open(tempmasking_file, 'rb') as pklfile:
    temporal_masking_fn = dill.load(pklfile)

spatial_unmasking_fn = pd.read_csv(spatial_unmasking_file)

kwargs['temporal_masking_thresholds'] = temporal_masking_fn
kwargs['spatial_release_fn'] = np.array(spatial_unmasking_fn)[:,1:]
   
kwargs['call_directionality'] = lambda  X, A=7 : A*(np.cos(np.deg2rad(X))-1)

kwargs['hearing_directionality'] = lambda X, B=2 : B*(np.cos(np.deg2rad(X))-1)


common_parameters = '../common_simulation_parameters.paramset'
with open(common_parameters,'wb') as paramsfile:
    dill.dump(kwargs, paramsfile)
