# -*- coding: utf-8 -*-
""" Wrapper that sets up a CPN simulation with command line interface

Created on Fri Aug 23 17:00:56 2019

@author: tbeleyur
"""
import argparse
import hashlib
import multiprocessing as mp
import dill 
from joblib import Parallel, delayed
import glob
import pdb
import uuid
import os 
import sys 
wd="../"
sys.path.append(wd+"/CPN/")
import numpy as np 

from the_cocktail_party_nightmare import run_CPN, No_situation_for_given_nd

def make_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

           
argparser = argparse.ArgumentParser()
argparser.add_argument('-param_file', action='store',dest='parameter_fileformat')
argparser.add_argument('-numCPUS','-numCPUs', action='store', dest='num_cpus',
                                   type=int, 
                                   default=mp.cpu_count())
argparser.add_argument('-name', action='store',
                       dest='name')
argparser.add_argument('-info', action='store',
                       dest='info', default='')
argparser.add_argument('-dest_folder', action='store',
                       dest='dest_folder', default='.')
args = argparser.parse_args()

def load_parameters(path_to_parameters):
    '''
    '''
    try:
        with open(path_to_parameters, 'rb') as paramsfile:
            parameter_instance = dill.load(paramsfile)
        return(parameter_instance)
    except:
        raise


def load_all_parameters(multi_parameter_files):
    '''
    '''
    all_parameter_sets = map(load_parameters, multi_parameter_files)
    return(all_parameter_sets)

def save_simulation_outputs(dest_folder,
                            simulation_identifiers, sim_output):
    '''
    Parameters
    ----------
    
    dest_folder : 
    simulation_identifiers : list with information required to identify 
                             and re-run a simulation if necessary.    

    sim_output : list with 3 objects inside
                echo_ids 
                sounds_in_ipi
                group_geometry 
    '''
    picklefilename = '_'.join([simulation_identifiers['name'],
                              simulation_identifiers['uuid'],
                              str(simulation_identifiers['np.random.seed']),
                              '.simresults'])
    full_filepath = os.path.join(dest_folder, picklefilename)
    try:
        with open(full_filepath, 'wb') as picklefile:
            dill.dump([simulation_identifiers, sim_output], picklefile)
        return(True)
    except:
        raise IOError(picklefilename + ' not saved!!')

def run_and_save_one_simulation(id_and_params):
    '''
    Runs one simulation and saves the output as a Pickle file. 

    Parameters
    ----------
    
    id_and_params : list with two entries :
        simulation_identifiers : dictionary with information to uniquely identify
                             each and every simulation run.

        parameter_set : dictionary with parameters to initialise a simulation

    Returns
    --------
    success : Boolean. 
              True if the run + saving was succesful , False if not. 
    
    Note
    -------
    
    
    '''
    simulation_identifiers, parameter_set = id_and_params
    #generate unique id for this simulation 
    file_uuid = str(uuid.uuid4())
    unique_seed = int(hashlib.sha1(file_uuid.encode('utf-8')).hexdigest(), 16) % (2**30)
    np.random.seed(unique_seed)
    # run simulation 
    try:
        num_echoes, sim_output = run_CPN(**parameter_set)
    
    # print(sim_output)
    # save outputs 
        simulation_identifiers['uuid'] = file_uuid
        simulation_identifiers['np.random.seed'] = unique_seed

        success = save_simulation_outputs(parameter_set['dest_folder'],
                                            simulation_identifiers, sim_output)
    except No_situation_for_given_nd as error:
        print(str(error)) # Very bad mistake
        success=False
    if not success:
        print('Simulation ' , file_uuid, 'could not be saved')
    return(success)

def run_multiple_simulations(name, info,
                                 parameter_file_format,
                             num_CPUs,
                             dest_folder, varying_param):
    ''' Set up the simulation to run parallely
    
    Parameters
    ----------
    name : str. 
           A *short* descriptive name to identify what this set of simulation runs 
           was doing. 

    info : str. 
           A somewhat longer description of what this set of simulation runs was
           doing. 

    Nruns : int>0. 
            Number of simulation repicates to run

    parameter_file_format : string which identifies the parameter file formats
                        with parameters to initialise and run 
                        the simulations
    num_CPUs: int. 
              Number of CPUs to run the code on. 
              

    dest_folder : str or path.
                  where the simulation outputs will be saved.   
    Returns
    --------
    success : list. 
              Nrun Boolean entries. True if the run+ saving was succesful
              False if not. 
            
    '''

    parameter_files = glob.glob(parameter_file_format)
    # print(parameter_files)
    parameter_sets = load_all_parameters(parameter_files)
    
    id_and_paramsets = []
    for i,each_paramset in enumerate(parameter_sets):
        folder_id= str(each_paramset[varying_param])
        simulation_identifiers = {}
        simulation_identifiers['name'] = name
        simulation_identifiers['info'] = info
        simulation_identifiers['parameter_set'] = each_paramset
        make_dir(dest_folder+folder_id)
        each_paramset['dest_folder'] = dest_folder + folder_id
        for _ in range(each_paramset['Nruns']):
            id_and_paramsets.append([simulation_identifiers, each_paramset])
        print(dest_folder+folder_id)
    
    success = Parallel(n_jobs=num_CPUs,
                           verbose=1, backend="loky")(map(delayed(run_and_save_one_simulation),
                                                     id_and_paramsets))
    return(success)

    
class FailedParameterLoading(ValueError):
    pass

# if __name__  == '__main__':
#     run_multiple_simulations(args.name, args.info, 
#                              args.parameter_fileformat,
#                              args.num_cpus, args.dest_folder)
# print(os.getcwd())
# parameter_file_format="./effect_of_heading_direction/paramsets/*"
# print(glob.glob(parameter_file_format)[0][86:][:-13])

# run_multiple_simulations(name="heading_direction", info="trail run",
#                                 parameter_file_format="./effect_of_heading_direction/paramsets/*",
#                             num_CPUs=10,
#                             dest_folder="./effect_of_heading_direction/store_results/", varying_param= "heading_variation")       
    
# run_multiple_simulations(name="focal_bat_r", info="trail run",
#                                 parameter_file_format="./effect_of_position_focal_bat_r/paramsets/*",
#                             num_CPUs=10,
#                             dest_folder="./effect_of_position_focal_bat_r/store_results/", varying_param="central_bat")       

# run_multiple_simulations(name="focal_bat_theta", info="trail run",
#                                 parameter_file_format="./effect_of_position_focal_bat_theta/paramsets/*",
#                             num_CPUs=10,
#                             dest_folder="./effect_of_position_focal_bat_theta/store_results/", varying_param="noncentral_bat")       

run_multiple_simulations(name="focal_bat_d", info="trail run",
                                parameter_file_format="./effect_of_position_focal_bat_d/paramsets/*",
                            num_CPUs=20,
                            dest_folder="./effect_of_position_focal_bat_d/store_results/", varying_param= "change_focal_bat_d")      