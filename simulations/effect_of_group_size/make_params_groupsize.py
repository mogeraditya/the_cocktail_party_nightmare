#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Create simulation parameters and run to study the effect of group size 

Created on Sat Aug 24 22:10:10 2019

@author: tbeleyur
"""

import glob
import os
import sys
import dill 
import pandas  as pd
import numpy as np 
import statsmodels.api as sm


join_into_string = lambda Y: '*'.join(map(lambda X:str(X), Y))

wd="/home/adityamoger/Documents/GitHub/cocktail_clone/"

all_group_sizes = [5,10,30,50,75,100,200]
interpulse_duration = [0.1]
call_duration = [0.0025]
shadowing = [True]
source_level = [100]
spacing = [0.5]
group_heading_variation = [10]
atmospheric_absorption = [-1]
number_of_simulation_runs =  1

with open(wd+'/common_simulation_parameters.paramset','rb') as pklfile:
    simulation_parameters = dill.load(pklfile)
    
simulation_parameters['echoes_beyond_ipi'] = True

i = 0
for group_size in all_group_sizes:
    for ipi in interpulse_duration:
        for call_durn in call_duration:
            for w_shadowing in shadowing:
                for SL in source_level:
                    for interbat_spacing in spacing:
                        for heading_variation in group_heading_variation:
                            for atm_abs in atmospheric_absorption:
                                i += 1
                                simulation_parameters['Nbats']  = group_size
                                simulation_parameters['Nruns']  = number_of_simulation_runs
                                description = str(simulation_parameters['Nruns']) + 'runs across'+str(group_size)
                                simulation_parameters['detailed description'] = description
                                simulation_parameters['interpulse_interval'] = ipi
                                simulation_parameters['echocall_duration'] = call_durn
                                simulation_parameters['implement_shadowing'] = w_shadowing
                                simulation_parameters['source_level'] = {'dBSPL' : SL, 
                                                                                   'ref_distance':1.0}
                                simulation_parameters['min_spacing'] = interbat_spacing
                                simulation_parameters['heading_variation'] = heading_variation
                                simulation_parameters['atmospheric_attenuation'] = atm_abs
    
                                all_params = [group_size, 
                                              ipi,call_durn,w_shadowing,SL,interbat_spacing,heading_variation,atm_abs]
                                variables_as_string = join_into_string(all_params)
    
                                param_filename = 'simulation_parameters_' + variables_as_string+'_.paramset'
                                os.chdir("/home/adityamoger/Documents/GitHub/cocktail_clone/simulations/effect_of_group_size/paramsets/central/")
                                with open(param_filename,'wb') as pklfile:
                                    dill.dump(simulation_parameters, pklfile)
