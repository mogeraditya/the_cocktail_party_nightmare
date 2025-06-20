
import glob
import os
import sys

import dill 
import pandas  as pd
import numpy as np 
import statsmodels.api as sm


join_into_string = lambda Y: '*'.join(map(lambda X:str(X), Y))

group_sizes = [5,10,30,50]
interpulse_duration =  0.1
call_duration =  0.0025
shadowing = True
source_level =  100
spacing = 0.5
group_heading_variation = [0, 10, 30, 50, "swarming", "flocking"]
atmospheric_absorption =  -1
number_of_simulation_runs = 50

# print(os.getcwd())
# os.chdir("../../")
# print(os.getcwd())
with open('../../common_simulation_parameters.paramset','rb') as pklfile:
    simulation_parameters = dill.load(pklfile)
    
simulation_parameters['echoes_beyond_ipi'] = True

i = 0

for group_size in group_sizes:
    for heading_var in group_heading_variation:
        i += 1; 
        simulation_parameters['Nbats']  = group_size
        simulation_parameters['Nruns']  = number_of_simulation_runs
        description = str(simulation_parameters['Nruns']) + 'runs across'+str(group_size)
        simulation_parameters['detailed description'] = description
        simulation_parameters['interpulse_interval'] = interpulse_duration
        simulation_parameters['echocall_duration'] = call_duration
        simulation_parameters['implement_shadowing'] = shadowing
        simulation_parameters['source_level'] = {'dBSPL' : source_level, 
                                                           'ref_distance':1.0}
        simulation_parameters['min_spacing'] = spacing
        simulation_parameters['heading_variation'] = heading_var
        simulation_parameters['atmospheric_attenuation'] = atmospheric_absorption
        
        # simulation_parameters['noncentral_bat'] = focal_bat_position

        all_params = [interpulse_duration, call_duration, shadowing,
                          source_level, spacing,
                          atmospheric_absorption, group_size,
                          heading_var]
        variables_as_string = join_into_string(all_params)

        param_filename = 'simulation_parameters_' + variables_as_string+'_.paramset'

        with open("../effect_of_heading_direction/paramsets/"+param_filename,'wb') as pklfile:
            dill.dump(simulation_parameters, pklfile)
