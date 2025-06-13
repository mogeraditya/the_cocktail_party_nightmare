import glob
import os
import sys

import dill 
import pandas  as pd
import numpy as np 
import statsmodels.api as sm
import math
import matplotlib.pyplot as plt

def make_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

wd="/home/adityamoger/Documents/GitHub/cocktail_clone/"
sys.path.append(wd+"/CPN/")
import the_cocktail_party_nightmare as CPN
import create_flocking_swarming_angles as flsw


join_into_string = lambda Y: '*'.join(map(lambda X:str(X), Y))

group_sizes = [75]
interpulse_duration =  0.1
call_duration =  0.0025
shadowing = True
source_level =  100
spacing = 0.5
group_heading_variation = 10
atmospheric_absorption =  -1
number_of_simulation_runs =  200
radial_distance = [1]
azimuth_location = [np.pi/4, 5*np.pi/4] 

with open(wd+'/common_simulation_parameters.paramset','rb') as pklfile:
    simulation_parameters = dill.load(pklfile)
    
simulation_parameters['echoes_beyond_ipi'] = True


def compute_u_v_for_plotting(points, angles):
    u_v_s= []; i=0
    for angle in angles:
        point= points[i]
        angle_in_radians= math.radians(angle)
        U= np.cos(angle_in_radians)
        V= np.sin(angle_in_radians)
        
        u_v_s.append([point[0], point[1], U, V])
        i+=1

    return np.array(u_v_s)

i = 0
for group_size in group_sizes:
    for r in radial_distance:
        for theta in azimuth_location:
            i += 1
            focal_bat_position = (r,theta)
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
            simulation_parameters['noncentral_bat'] = focal_bat_position

            bat_positions, bats_orientations = CPN.place_bats_inspace(**simulation_parameters)
            # nearby, focal = CPN.generate_surroundpoints_w_poissondisksampling(npoints=simulation_parameters['Nbats'],
            #                                                        nbr_distance= simulation_parameters['min_spacing'],
            #                                                        **simulation_parameters)


            all_bat_positions= np.row_stack((bat_positions[1], bat_positions[0]))
            flock= flsw.flocking(all_bat_positions, 32)
            u_v_s= compute_u_v_for_plotting(all_bat_positions, flock)
            for vector in u_v_s:
                plt.quiver(*vector, angles="uv")
            plt.plot(all_bat_positions[1:][:,0], all_bat_positions[1:][:,1], ".")
            plt.plot(all_bat_positions[0][0], all_bat_positions[1][1], ".", color="red")
            plt.show()

            swarm= flsw.swarming(**simulation_parameters)
            u_v_s= compute_u_v_for_plotting(all_bat_positions, swarm)
            for vector in u_v_s:
                plt.quiver(*vector, angles="uv")
            plt.plot(all_bat_positions[1:][:,0], all_bat_positions[1:][:,1], ".")
            plt.plot(all_bat_positions[0][0], all_bat_positions[1][1], ".", color="red")
            plt.show()





