import glob
import os
import sys

import dill 
import pandas  as pd
import numpy as np 
import statsmodels.api as sm

def make_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

wd="../../"

join_into_string = lambda Y: '*'.join(map(lambda X:str(X), Y))

group_sizes = [100,75,50,40,30,20,10,5]
interpulse_duration =  0.1
call_duration =  0.0025
shadowing = True
source_level =  100
spacing = 0.5
group_heading_variation = 10
atmospheric_absorption =  -1
number_of_simulation_runs = 50
# radial_distance = [0.2, 0.5, 0.75]
azimuth_location = [1] 
radial_values=[np.pi/2]
threshold= 0.25
with open(wd+'/common_simulation_parameters.paramset','rb') as pklfile:
    simulation_parameters = dill.load(pklfile)
    
simulation_parameters['echoes_beyond_ipi'] = True

i = 0

for group_size in group_sizes:
    for radial_value in radial_values:
    # for r in radial_distance:
        for theta in azimuth_location:
            i += 1
            # focal_bat_position = (0,0)
            simulation_parameters["change_focal_bat_d"] = (radial_value, theta, threshold)
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
            simulation_parameters['heading_variation'] = group_heading_variation
            simulation_parameters['atmospheric_attenuation'] = atmospheric_absorption
            # simulation_parameters['central_bat'] = focal_bat_position

            all_params = [group_size,interpulse_duration, call_duration, shadowing,
                            source_level, spacing,
                            group_heading_variation,
                            atmospheric_absorption, group_size,
                            theta, threshold, radial_value]
            variables_as_string = join_into_string(all_params)

            param_filename = 'simulation_parameters_' + variables_as_string+'_.paramset'
            make_dir("./paramsets/")
            with open("./paramsets/"+param_filename,'wb') as pklfile:
                dill.dump(simulation_parameters, pklfile)
