import random
import dill 
import datetime as dt
import glob 
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import pandas as pd
import scipy.spatial as spatial
import sys 
wd=".."
sys.path.append(wd+'//CPN//')
import numpy as np 
import statsmodels.api as sm
from tqdm import tqdm_notebook, tqdm
import brokenaxes
import math
import os
from analysis_misc_functions import *
import IPython

group_heading_variation = [0, 10, 30, 50, "swarming", "milling"]

for heading_var in group_heading_variation:
    dir_to_store_plots= "./plots_heading_radial/"+str(heading_var)
    make_dir(dir_to_store_plots)
    results_folder = wd+"/simulations/effect_of_heading_plus_radial/store_results/swarming/"
    list_of_folder_ids= glob.glob(results_folder+"*")
    list_of_folder_ids= [i[len(results_folder):] for i in list_of_folder_ids]
    list_of_folder_ids= np.sort(list_of_folder_ids)
    list_of_folder_ids_pruned= []
    for folder_id in list_of_folder_ids:
        all_results = glob.glob(results_folder+str(folder_id)+'/*.simresults')
        if len(all_results)==0:
            continue
        list_of_folder_ids_pruned.append(folder_id)

    list_of_folder_ids_array_version= [i.split(",") for i in list_of_folder_ids_pruned]
    list_of_folder_ids_array_version= [(n[0][1:],n[1][1:]) for n in list_of_folder_ids_array_version]
    list_of_folder_ids_array_version= np.array([(float(n[0]), float(n[1])) for n in list_of_folder_ids_array_version])

    group_sizes= [5,10,20,30,50,75]
    list_of_folder_ids_array_version
    folder_id= list_of_folder_ids[0]
    all_results = glob.glob(results_folder+folder_id+'/*.simresults')
    some_results =random.sample(all_results, int(len(all_results)*1.0))
    extraction_fns = [get_num_echoes_heard, get_group_size, get_furthest_bat2bat_distance,
                    get_nearest_neighbour_distances, get_detection_distance, get_run_uuid,
                    get_run_random_seed,make_paramset_id, extract_parameter_values,
                    get_echo_levels, get_detection_azimuth]
    keyword_arguments = {'nearest_nbrs':3}
    keyword_arguments['variables_to_extract'] = ['heading_variation', 'echocall_duration','atmospheric_attenuation','min_spacing','source_level',
                            'interpulse_interval', 'implement_shadowing',
                                                ]
    what_to_do = input('Do you want to "rerun" from scratch or 2, "reanalyse" from previous data')
    task_type= input('What type of analysis is this?')
    what_to_do
    if what_to_do == 'rerun':
        array_containing_all_simulation_datasets_raw = []
        iterant=0
        yyyymmdd = dt.datetime.now()
        timestamp = str([yyyymmdd.year,yyyymmdd.month,yyyymmdd.day,yyyymmdd.hour])
        make_dir(dir_to_store_plots+"/"+task_type+'_data'+timestamp+"/")
        for folder_id in list_of_folder_ids_pruned:
            all_results = glob.glob(results_folder+str(folder_id)+'/*.simresults')
            some_results =random.sample(all_results, int(len(all_results)*1.0))
            if len(some_results)==0:
                continue
            
            %time extracted_simdata = Parallel(n_jobs=4)(delayed(load_and_extract)(each, extraction_fns, **keyword_arguments) for each in tqdm(some_results))
            # print(extracted_simdata)
            echoes_heard = []
            group_size = []
            group_diameter = []
            nearest_3nbrs = []
            nbr_detection_range = []
            uuid = []
            seed = []
            param_ids = []
            parameter_values = []
            echo_levels = []
            detection_angle = []
            for each in extracted_simdata:
                for variable, list_to_append in zip(each, [echoes_heard, group_size, group_diameter,
                                                        nearest_3nbrs, nbr_detection_range,
                                                        uuid, seed, param_ids, parameter_values,
                                                        echo_levels, detection_angle]):
                    list_to_append.append(variable)
            simulation_data = pd.DataFrame(data = {'nbrs_detected':echoes_heard,
                                            'group_size':group_size,
                                            'group_diameter':group_diameter,
                                            'nearest_neighbour_distance':nearest_3nbrs,
                                            'nbrs_detected_distance':nbr_detection_range,
                                            'uuid':uuid,
                                            'seed':seed,
                                            'paramset_id':param_ids,
                                            'parameters_joint':parameter_values,
                                            'echo_levels':echo_levels,
                                            'detection_azimuth':detection_angle})
            column_names = keyword_arguments['variables_to_extract']

            for each in column_names:
                simulation_data[each] = np.nan

            %time simulation_data = simulation_data.apply(split_parameters_joint_to_separate_columns, 1, col_names=column_names)    
            

            simulation_data["focal_bat"]= len(simulation_data)*[list_of_folder_ids_array_version[iterant]]

            ##### thanks to https://realpython.com/fast-flexible-pandas/#but-i-heard-that-pandas-is-slow
            data_store = pd.HDFStore(dir_to_store_plots+"/"+task_type+'_data'+timestamp+'/'+folder_id+'.h5')
            data_store['simulation_data'] = simulation_data
            data_store.close()

            array_containing_all_simulation_datasets_raw.append(simulation_data)
            iterant+=1
            
    if what_to_do == 'reanalyse':
        # show all datasets in the folder
        data_stores_in_folder = glob.glob(task_type+'_*]')

        datastores_by_index = pd.DataFrame(data={'data_stores':data_stores_in_folder})
        print(datastores_by_index)
        # user-inputs dataset choice:
        which_dataset = int(input('Give the index of the dataset you want to'))
        print()
        # load a pre-exisiting dataset
        #load the saved data 
        list_of_h5_files_in_folder= making_glob_work_properly(datastores_by_index["data_stores"][which_dataset])
        # print(list_of_h5_files_in_folder); print("./"+datastores_by_index["data_stores"][which_dataset]+'/*')
        array_containing_all_simulation_datasets_raw = []
        list_of_folder_ids_pruned, list_of_folder_ids_array_version= run_while_reanalysis_to_change_folder_ids(list_of_h5_files_in_folder)
        for h5_file in list_of_h5_files_in_folder:
            data_load = pd.HDFStore(dir_to_store_plots+"/"+datastores_by_index["data_stores"][which_dataset]+'/'+h5_file)
            simulation_data = data_load['simulation_data']
            array_containing_all_simulation_datasets_raw.append(simulation_data)
            data_load.close()

    # else:
    #     raise IOError('Invalid input given!!')
    array_containing_all_simulation_datasets_pruned= []
    for group_size in group_sizes:
        for simulation_data in array_containing_all_simulation_datasets_raw:
            indices_given_group_size= simulation_data["group_size"]==group_size
            if len(simulation_data[indices_given_group_size])<100:
                array_containing_all_simulation_datasets_pruned.append(simulation_data[~indices_given_group_size])

    array_containing_all_simulation_datasets=[]

    for simulation_data in array_containing_all_simulation_datasets_raw:
        if False: #simulation_data["focal_bat"][0][0]==5.0 or simulation_data["focal_bat"][0][0]==6.0:
            continue
        else:
            array_containing_all_simulation_datasets.append(simulation_data)



    unique_angles_in_dataset= np.unique(list_of_folder_ids_array_version[:,1])
    list_containing_simulated_data_sorted_by_angles= {}
    for angle in unique_angles_in_dataset:
        list_containing_simulated_data_sorted_by_angles[angle]= []
        sublist_true_false_labels= list_of_folder_ids_array_version[:,1]==angle
        sublist_given_angle= passing_true_false_indices_through_lists(array_containing_all_simulation_datasets_raw, sublist_true_false_labels)
        list_containing_simulated_data_sorted_by_angles[angle].append(sublist_given_angle)
    # %matplotlib inline
    i=0
    for angle in list_containing_simulated_data_sorted_by_angles.keys():
        print(angle)
        subset_given_angle= list_containing_simulated_data_sorted_by_angles[angle][0]

        for simulation_data in subset_given_angle:
            if len(simulation_data)!=0:
                plt.plot(simulation_data['group_size'],
                            simulation_data['group_diameter'], 
                        '*',label=simulation_data["focal_bat"][0],alpha=0.2)
            i+=1
        plt.legend()
        # make_dir(dir_to_store_plots+"/test_accuracy_of_data.png")
        plt.savefig(dir_to_store_plots+"/test_accuracy_of_data.png")
        plt.clf()
    simulation_data.keys()
    by_minspacing = simulation_data.groupby('min_spacing')
    nearest_nbr_dist = {}
    for spacing, df in by_minspacing:
        nearest_nbr_distances = df['nearest_neighbour_distance'].reset_index(drop=True)
        nearest_nbr_dist[spacing] = np.concatenate(nearest_nbr_distances)
    nearest_nbr_dists_for_plot = [nearest_nbr_dist[0.5]]
    iterant1_axis=0
    # fig, axs = plt.subplots(1,3, figsize=(27,9))


    for angle in list_containing_simulated_data_sorted_by_angles.keys():
            subset_given_angle= list_containing_simulated_data_sorted_by_angles[angle][0]; iterant2_color=0
            colors=["blue", "red", "green", "black"]
            plt.figure(figsize=(21,10))
            i=1
            for simulation_data in subset_given_angle:
                    data_by_groupsize = simulation_data.groupby(['group_size'])
                    nbrs_detected_data = {}
                    nbrs_detected_data['num_nbrs'] = []
                    nbrs_detected_data['median_num_nbrs'] = []
                    nbrs_detected_data['90%ile_num_nbrs'] = []
                    groupsizes= []
                    for groupsize, dataframe in data_by_groupsize:
                            groupsizes.append(groupsize)
                            nbrs_detected_data['num_nbrs'].append(np.array(dataframe['nbrs_detected']))
                            nbrs_detected_data['median_num_nbrs'].append(np.median(np.array(dataframe['nbrs_detected'])))
                            nbrs_detected_data['90%ile_num_nbrs'].append(np.percentile(np.array(dataframe['nbrs_detected']), 90))

                    breakup_x_axis  = np.arange(1,len(groupsizes)+1)+0.5
                    nbrs_detected_data['median_num_nbrs']

                    # plt.figure(figsize=(7,3.42))
                    plt.subplot(2,3,i)
                    plt.violinplot(nbrs_detected_data['num_nbrs'][:], showextrema=False)
                    plt.plot(np.arange(1,len(groupsizes)+1),nbrs_detected_data['median_num_nbrs'][:],
                            '--*',linewidth=2, alpha=0.5,label="median, r="+str(simulation_data["focal_bat"][0]))#'median')
                    # print(np.arange(1,len(groupsizes)))
                    plt.plot(np.arange(1,len(groupsizes)+1), nbrs_detected_data['90%ile_num_nbrs'][:],
                            '-..',linewidth=2, alpha=0.5, label="90% percentile, r="+str(simulation_data["focal_bat"][0]))#label='90%ile')
                    # plt.vlines(breakup_x_axis, -1,20, 'w',linewidth=30,alpha=1.0)
                    plt.yticks(np.arange(0,21,1))
                    plt.xticks(range(1,len(groupsizes)+1),groupsizes[:])
                    plt.legend(fontsize=8)
                    plt.ylabel('Neighbours detected \n per call emission', fontsize=10)
                    plt.xlabel('Group size', fontsize=10)
                    # # plt.ylim(0,20)

                    # axs[i].tight_layout()
                    plt.title("neighbours detected when focal bat is at angle "+str(np.degrees(np.pi/2)))
                    iterant2_color+=1; i+=1
            make_dir(dir_to_store_plots+f"/no_of_neighbours/")
            plt.savefig(dir_to_store_plots+f"/no_of_neighbours/{angle}.png")
            plt.clf()        
            iterant1_axis+=1
                    # plt.savefig('Figure2_groupsize_1200dpi.pdf',
                    #         bbox_inches='tight', pad_inches=0.05, format='pdf', dpi=1200)
            
            # print(groupsizes); print(nbrs_detected_data['median_num_nbrs'])
    
    
    
    # Verticle transect

    # moving from $(0.75, \frac {3pi}{2})$ to $(0.75, \frac {pi}{2})$ [r, theta] is similar to moving in a ($x=0$) transect centered at the centermost bat.
    # Note: All the heading vectors are in 90 degree direction ( $\pm$ 10 degree)
    transect_subsets= [[6, 3*np.pi/2],[5, 3*np.pi/2],[4, 3*np.pi/2],[3, 3*np.pi/2],[2, 3*np.pi/2],[1, 3*np.pi/2],[1, np.pi/2],[2, np.pi/2], [3, np.pi/2], [4, np.pi/2], [5, np.pi/2], [6, np.pi/2]]
    iterant_transect=0
    sublist_of_simulation_data_intransect=[]
    for r_theta in transect_subsets:
        for simulation_data in array_containing_all_simulation_datasets:
            if (simulation_data["focal_bat"][0]==r_theta).all():
                sublist_of_simulation_data_intransect.append(simulation_data)    

    sign_based_on_distance_from_centermost_point= [-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1]     

    len(sublist_of_simulation_data_intransect)
    # len(transect_subsets)
    colors=["blue", "red", "green", "black", "orange", "cyan"]

    fig, axs = plt.subplots(1,1, figsize=(27,9))
    iterant2_color=0
    for simulation_data in sublist_of_simulation_data_intransect:
            data_by_groupsize = simulation_data.groupby(['group_size'])
            nbrs_detected_data = {}
            nbrs_detected_data['num_nbrs'] = []
            nbrs_detected_data['median_num_nbrs'] = []
            nbrs_detected_data['90%ile_num_nbrs'] = []
            groupsizes= []
            for groupsize, dataframe in data_by_groupsize:
                    groupsizes.append(groupsize)
                    nbrs_detected_data['num_nbrs'].append(np.array(dataframe['nbrs_detected']))
                    nbrs_detected_data['median_num_nbrs'].append(np.median(np.array(dataframe['nbrs_detected'])))
                    nbrs_detected_data['90%ile_num_nbrs'].append(np.percentile(np.array(dataframe['nbrs_detected']), 90))

            breakup_x_axis  = np.arange(1,len(groupsizes)+1)+0.5
            nbrs_detected_data['median_num_nbrs']

            # plt.figure(figsize=(7,3.42))
            # plt.subplot(1,2,i)
            
            axs.violinplot(nbrs_detected_data['num_nbrs'][:], showextrema=False)
            # axs.plot(np.arange(1,len(groupsizes)+1),nbrs_detected_data['median_num_nbrs'][:],
            #         '--*',linewidth=2, alpha=0.5,label="median, r="+str(sign_based_on_distance_from_centermost_point[iterant2_color]*simulation_data["focal_bat"][0][0]), 
            #         color= colors[iterant2_color])#'median')
            # print(np.arange(1,len(groupsizes)))
            axs.plot(np.arange(1,len(groupsizes)+1), nbrs_detected_data['90%ile_num_nbrs'][:],
                    '-..',linewidth=2, alpha=0.5, label="90% percentile, r="+str(sign_based_on_distance_from_centermost_point[iterant2_color]*simulation_data["focal_bat"][0][0])
                    )#label='90%ile')
            # plt.vlines(breakup_x_axis, -1,20, 'w',linewidth=30,alpha=1.0)
            axs.set_yticks(np.arange(0,21,1))
            axs.set_xticks(range(1,len(groupsizes)+1),groupsizes[:])
            axs.legend(fontsize=12)
            axs.set_ylabel('Neighbours detected \n per call emission', fontsize=10)
            axs.set_xlabel('Group size', fontsize=10)
            # # plt.ylim(0,20)

            # axs[i].tight_layout()
            axs.set_title("neighbours detected when focal bat is at angle "+str(np.degrees(simulation_data["focal_bat"][0][1])))
            iterant2_color+=1
    plt.show()
    plt.clf()

    # separate the violin plots out for different group sizes to see clear differences between distance from center

    colors=["blue", "red", "green", "black", "pink", "cyan", "yellow"]

    # fig, axs = plt.subplots(1,1, figsize=(27,9))
    i=0; 
    number_of_columns=3; number_of_rows= int(np.ceil(len(group_sizes)/number_of_columns))
    plt.figure(figsize=(20,12))
    for group_size in group_sizes:
        nbrs_detected_data = {}
        nbrs_detected_data['num_nbrs'] = []
        nbrs_detected_data['median_num_nbrs'] = []
        nbrs_detected_data['90%ile_num_nbrs'] = []
        
        x_ticks= []; iterant_sign=-1
        for simulation_data in sublist_of_simulation_data_intransect:
            iterant_sign+=1

            try:
                dataframe = simulation_data.groupby(['group_size']).get_group(group_size)    
                if len(dataframe["nbrs_detected"])!=100:
                    continue
                nbrs_detected_data['num_nbrs'].append(np.array(dataframe['nbrs_detected']))
                nbrs_detected_data['median_num_nbrs'].append(np.median(np.array(dataframe['nbrs_detected'])))
                nbrs_detected_data['90%ile_num_nbrs'].append(np.percentile(np.array(dataframe['nbrs_detected']), 90))    
            except KeyError:
                continue
            x_ticks.append(sign_based_on_distance_from_centermost_point[iterant_sign]*simulation_data["focal_bat"][0][0])
        

        plt.subplot(number_of_rows,number_of_columns,i+1)
        plt.violinplot(nbrs_detected_data['num_nbrs'], showmeans=True)
        # plt.boxplot(nbrs_detected_data['num_nbrs'], showmeans=True)
        plt.title(f"neighbours detected when focal bat is at angle {90}; groupsize= {group_size}", fontsize= "10")
        plt.xticks(np.arange(1,len(x_ticks)+1),x_ticks)
        # plt.show()
        plt.ylim(0,group_size+0.5)
        plt.axhline(group_size-1, ls="-.", color="red", label= "group_size - 1")
        plt.xlabel("y_distance from centermostpoint", fontsize= "8")
        plt.ylabel("number of neighburs detected", fontsize= "8")
        plt.legend(fontsize= "8")

        i+=1
    label= "/number_of_detected_neighbours/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_groupsize.png")
    plt.clf() 
    # Zooming into the plots to see finer differences between the diffrent violin plots


    # fig, axs = plt.subplots(1,1, figsize=(27,9))
    i=0; 
    number_of_columns=3; number_of_rows= int(np.ceil(len(group_sizes)/number_of_columns))
    plt.figure(figsize=(20,12))
    for group_size in group_sizes:
        nbrs_detected_data = {}
        nbrs_detected_data['num_nbrs'] = []
        nbrs_detected_data['median_num_nbrs'] = []
        nbrs_detected_data['90%ile_num_nbrs'] = []
        
        x_ticks= []; iterant_sign=-1
        for simulation_data in sublist_of_simulation_data_intransect:
            iterant_sign+=1
            try:
                dataframe = simulation_data.groupby(['group_size']).get_group(group_size)    
                if len(dataframe["nbrs_detected"])<90:
                    continue
                nbrs_detected_data['num_nbrs'].append(np.array(dataframe['nbrs_detected']))
                nbrs_detected_data['median_num_nbrs'].append(np.median(np.array(dataframe['nbrs_detected'])))
                nbrs_detected_data['90%ile_num_nbrs'].append(np.percentile(np.array(dataframe['nbrs_detected']), 90))    
            except KeyError:
                continue
            x_ticks.append(sign_based_on_distance_from_centermost_point[iterant_sign]*simulation_data["focal_bat"][0][0])

            

        plt.subplot(number_of_rows,number_of_columns,i+1)
        plt.violinplot(nbrs_detected_data['num_nbrs'], showmeans=True)
        # plt.boxplot(nbrs_detected_data['num_nbrs'], showmeans=True)
        plt.title(f"neighbours detected when focal bat is at angle {90}; groupsize= {group_size}", fontsize= "10")
        plt.xticks(np.arange(1,len(x_ticks)+1),x_ticks)
        # plt.show()
        # plt.ylim(0,group_size+0.5)
        # plt.axhline(group_size-1, ls="-.", color="red", label= "group_size - 1")
        plt.xlabel("y_distance from centermostpoint", fontsize= "8")
        plt.ylabel("number of neighburs detected", fontsize= "8")
        # plt.legend(fontsize= "8")

        i+=1
    label= "/number_of_detected_neighbours/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_radial.png")
    plt.clf() 
    # The number of neighbours detected shows an increase with group size initially and then drops sharply down beyond 30 bats. 

    ### How does the probability of detecting at least one neighbour change with group size?
    def calculate_geqNneighbour(num_detected_neighbours, N):
        '''Calculates the proportion of an array-like
        object that is >= N.
        
        '''
        geq_N = np.float64(np.sum(num_detected_neighbours>=N))
        if np.isnan(np.min(num_detected_neighbours)):
            return np.nan
        proportion_geq_N = geq_N/len(num_detected_neighbours)
        return(proportion_geq_N)
    simulation_data.keys()
    sublist_of_simulation_data_intransect[1]
    #store data
    store_probability_of_detection_n_neighbour= {}
    group_sizes= [5,10,20,30,50,75]
    i=0
    for simulation_data in sublist_of_simulation_data_intransect:
            data_by_groupsize = simulation_data.groupby(['group_size'])

            nbrs_detected_data = {}
            nbrs_detected_data["group_sizes"]= group_sizes
            nbrs_detected_data['num_nbrs'] = []
            nbrs_detected_data['median_num_nbrs'] = []
            nbrs_detected_data['90%ile_num_nbrs'] = []
            groupsizes= []

            for groupsize, dataframe in data_by_groupsize:
                    groupsizes.append(groupsize[0])
                    nbrs_detected_data['num_nbrs'].append(np.array(dataframe['nbrs_detected']))
                    nbrs_detected_data['median_num_nbrs'].append(np.median(np.array(dataframe['nbrs_detected'])))
                    nbrs_detected_data['90%ile_num_nbrs'].append(np.percentile(np.array(dataframe['nbrs_detected']), 90))
            nan_groupsizes= [i for i in group_sizes if i not in groupsizes]

            for nan_groupsize in nan_groupsizes:
                    nbrs_detected_data['num_nbrs']= [np.nan, *nbrs_detected_data['num_nbrs']]
                    nbrs_detected_data['median_num_nbrs']= [np.nan, *nbrs_detected_data['median_num_nbrs']]
                    nbrs_detected_data['90%ile_num_nbrs']= [np.nan, *nbrs_detected_data['90%ile_num_nbrs']]

            nbrs_detected_data['0neighbours'] = [ 1-calculate_geqNneighbour(each,1)  for each in nbrs_detected_data['num_nbrs']]
            nbrs_detected_data['geq_1neighbour'] = [ calculate_geqNneighbour(each,1)  for each in nbrs_detected_data['num_nbrs']]
            nbrs_detected_data['geq_2neighbour'] = [ calculate_geqNneighbour(each,2)  for each in nbrs_detected_data['num_nbrs']]
            nbrs_detected_data['geq_3neighbour'] = [ calculate_geqNneighbour(each,3)  for each in nbrs_detected_data['num_nbrs']]
            nbrs_detected_data['geq_4neighbour'] = [ calculate_geqNneighbour(each,4)  for each in nbrs_detected_data['num_nbrs']]

            nbrs_detected_data['detected_distances'] = [ np.concatenate(np.array(dataframe['nbrs_detected_distance'])) for groupsize, dataframe in data_by_groupsize]
            nbrs_detected_data['median_detected_distances'] = map(np.median, nbrs_detected_data['detected_distances'])
            nbrs_detected_data['90%ile_detected_distances'] = map(lambda X : np.percentile(X,90),
                                                            nbrs_detected_data['detected_distances'][:-1])

            nbrs_detected_data['neighbours_azimuth'] = [ np.concatenate(np.array(dataframe['detection_azimuth'])) for groupsize, dataframe in data_by_groupsize]
            nbrs_detected_data['median_neighbours_azimuth'] = [ np.median(each)  for each in nbrs_detected_data['neighbours_azimuth'] ]
            nbrs_detected_data['95%ile_neighbours_azimuth'] = [ np.percentile(each, 92.5)  for each in nbrs_detected_data['neighbours_azimuth'][:-1] ]
            nbrs_detected_data['5%ile_neighbours_azimuth'] = [ np.percentile(each, 2.5)  for each in nbrs_detected_data['neighbours_azimuth'][:-1] ]

            nbrs_detected_data['echolevels'] = [ np.concatenate(np.array(dataframe['echo_levels'])) for groupsize, dataframe in data_by_groupsize]
            nbrs_detected_data['median_echolevels'] = [ np.median(each) for each in nbrs_detected_data['echolevels']]
            nbrs_detected_data['5%ile_echolevels'] = [ np.percentile(each, 2.5) for each in nbrs_detected_data['echolevels'][:-1]]
            nbrs_detected_data['95%ile_echolevels'] = [ np.percentile(each, 97.5) for each in nbrs_detected_data['echolevels'][:-1]]

            for nan_groupsize in nan_groupsizes:
                    nbrs_detected_data['detected_distances']= [np.nan, *nbrs_detected_data['detected_distances']]
                    nbrs_detected_data['median_detected_distances']= [np.nan, *nbrs_detected_data['median_detected_distances']]
                    nbrs_detected_data['90%ile_detected_distances']= [np.nan, *nbrs_detected_data['90%ile_detected_distances']]

                    nbrs_detected_data['neighbours_azimuth']= [[np.nan], *nbrs_detected_data['neighbours_azimuth']]
                    nbrs_detected_data['median_neighbours_azimuth']= [np.nan, *nbrs_detected_data['median_neighbours_azimuth']]
                    nbrs_detected_data['5%ile_neighbours_azimuth']= [np.nan, *nbrs_detected_data['5%ile_neighbours_azimuth']]
                    nbrs_detected_data['95%ile_neighbours_azimuth']= [np.nan, *nbrs_detected_data['95%ile_neighbours_azimuth']]

                    nbrs_detected_data['echolevels']= [[np.nan], *nbrs_detected_data['echolevels']]
                    nbrs_detected_data['median_echolevels']= [np.nan, *nbrs_detected_data['median_echolevels']]
                    nbrs_detected_data['5%ile_echolevels']= [np.nan, *nbrs_detected_data['5%ile_echolevels']]
                    nbrs_detected_data['95%ile_echolevels']= [np.nan, *nbrs_detected_data['95%ile_echolevels']]
                    
            
            store_probability_of_detection_n_neighbour[simulation_data["focal_bat"][0][0]*sign_based_on_distance_from_centermost_point[i]]= nbrs_detected_data
            i+=1
    nan_groupsize
    i=0; 
    number_of_columns=3; number_of_rows= int(np.ceil(len(sublist_of_simulation_data_intransect)/number_of_columns))
    plt.figure(figsize= (21,21))
    for simulation_data in sublist_of_simulation_data_intransect:
            data_by_groupsize = simulation_data.groupby(['group_size'])
            groupsizes= group_sizes

            key_for_dict= simulation_data["focal_bat"][0][0]*sign_based_on_distance_from_centermost_point[i]
            nbrs_detected_data= store_probability_of_detection_n_neighbour[key_for_dict]

            plt.subplot(number_of_rows,number_of_columns,i+1)
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_1neighbour'], '-^', linewidth=6, alpha=0.5, label='$\geq1$ neighbour per call')
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_2neighbour'], '-s', linewidth=6, alpha=0.5, label='$\geq2$ neighbour per call')
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_3neighbour'], '-p', linewidth=6, alpha=0.5, label='$\geq3$ neighbour per call')
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_4neighbour'], '-H', linewidth=6, alpha=0.5, label='$\geq4$ neighbour per call')
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['0neighbours'], '-Dk', linewidth=6, alpha=0.5, label='No neighbours')

            plt.yticks(fontsize=10)
            plt.xticks(np.arange(1,len(groupsizes)+1),
                    groupsizes,fontsize=10)
            plt.legend()
            plt.ylabel('$P(detecting \geq \ X \ neighbours)$'+'\n'+ '$per \ call$'+str(simulation_data["focal_bat"][0][0]),
                    fontsize=10)
            plt.xlabel('Group size',
                    fontsize=10)
            plt.ylim(-0.05,1.05)
            
            i+=1
    label= "/probabilty_of_detection/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_groupsize.png")
    plt.clf() 
    nbrs_detected_data['geq_1neighbour']
    number_of_neighbours_atleast_detected= [0,1,2,3,4]
    i=0; 
    number_of_columns=1; number_of_rows= int(np.ceil(len(number_of_neighbours_atleast_detected)/number_of_columns))
    plt.figure(figsize= (21,50))

    for simulation_data in sublist_of_simulation_data_intransect:

            groupsizes= group_sizes

            key_for_dict= simulation_data["focal_bat"][0][0]*sign_based_on_distance_from_centermost_point[i]
            nbrs_detected_data= store_probability_of_detection_n_neighbour[key_for_dict]

            if sign_based_on_distance_from_centermost_point[i]<0:
                    ls= ".-"
            else:
                    ls= "--"

            plt.subplot(number_of_rows,number_of_columns,2)
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_1neighbour'], ls, alpha=0.5, label=f'distance={key_for_dict}', color=colors[int(simulation_data["focal_bat"][0][0]-1)])
            
            plt.subplot(number_of_rows,number_of_columns,3)
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_2neighbour'], ls, alpha=0.5, label=f'distance={key_for_dict}', color=colors[int(simulation_data["focal_bat"][0][0]-1)])
            
            plt.subplot(number_of_rows,number_of_columns,4)
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_3neighbour'], ls, alpha=0.5, label=f'distance={key_for_dict}', color=colors[int(simulation_data["focal_bat"][0][0]-1)])
            
            plt.subplot(number_of_rows,number_of_columns,5)
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['geq_4neighbour'], ls, alpha=0.5, label=f'distance={key_for_dict}', color=colors[int(simulation_data["focal_bat"][0][0]-1)])
            
            plt.subplot(number_of_rows,number_of_columns,1)
            plt.plot(np.arange(1,len(groupsizes)+1), 
                    nbrs_detected_data['0neighbours'], ls, alpha=0.5, label=f'distance={key_for_dict}', color=colors[int(simulation_data["focal_bat"][0][0]-1)])
            i+=1
    for i in range(1,6):
            if i==1:
                    plt.subplot(number_of_rows,number_of_columns,i)
                    plt.yticks(fontsize=10)
                    plt.xticks(np.arange(1,len(groupsizes)+1),
                            groupsizes,fontsize=10)
                    plt.legend()
                    plt.title(f"$\leq{i-1}$ neighbour per call;")
            else: 
                    plt.subplot(number_of_rows,number_of_columns,i)
                    plt.yticks(fontsize=10)
                    plt.xticks(np.arange(1,len(groupsizes)+1),
                            groupsizes,fontsize=10)
                    plt.legend()
                    plt.title(f"$\geq{i-1}$ neighbour per call;")
            plt.ylabel('$P(detecting \geq \ X \ neighbours)$'+'\n'+ '$per \ call$'+str(simulation_data["focal_bat"][0][0]),
                    fontsize=10)
            plt.xlabel('Group size',
                    fontsize=10)
            plt.ylim(-0.05,1.05)
    label= "/probabilty_of_detection/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_radial.png")
    plt.clf() 

    i=0; 
    number_of_columns=3; number_of_rows= int(np.ceil(len(sublist_of_simulation_data_intransect)/number_of_columns))
    plt.figure(figsize= (21,21))
    for simulation_data in sublist_of_simulation_data_intransect:
            groupsizes= group_sizes
            key_for_dict= simulation_data["focal_bat"][0][0]*sign_based_on_distance_from_centermost_point[i]
            nbrs_detected_data= store_probability_of_detection_n_neighbour[key_for_dict]

            plt.subplot(number_of_rows,number_of_columns,i+1)
            plt.boxplot(nbrs_detected_data['echolevels'][:], label=f"r={key_for_dict}");
            plt.vlines(range(1,len(nbrs_detected_data['95%ile_echolevels'])+1),nbrs_detected_data['5%ile_echolevels'],
                                                                            nbrs_detected_data['95%ile_echolevels'],
                                    label='95% data interval')

            plt.plot(range(1,len(groupsizes)+1), nbrs_detected_data['median_echolevels'], linewidth=6, alpha=0.5,
                    label='median detected echo level')
            plt.ylabel('Received level of detected echoes, \n dB SPL re 20$\mu Pa$',
                                                    fontsize=10)
            plt.xticks(np.arange(1,len(groupsizes)+1), groupsizes, fontsize=10)
            plt.hlines([20],1,len(groupsizes)+1)
            plt.yticks(np.arange(20,86,6),np.arange(20,86,6), fontsize=10)
            plt.text(1,21,'Hearing threshold', fontsize=10)
            plt.xlabel('Group size', fontsize=10)
            plt.legend()
            print(key_for_dict)
            i+=1
    label= "/received_level_dB/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_groupsize.png")
    plt.clf() 
    
    i=0; 
    d_distances_along_transect= [-4,-3,-2,-1,1,2,3,4]
    number_of_columns=3; number_of_rows= int(np.ceil(len(group_sizes)/number_of_columns))
    plt.figure(figsize= (21,21))
    for group_size in group_sizes:
            array_storing_info_for_plot= []
            iterant_sign=-1
            x_ticks=[]
            for key in d_distances_along_transect:
                    iterant_sign+=1
                    nbrs_detected_data= store_probability_of_detection_n_neighbour[key]
                    data_to_store= nbrs_detected_data['echolevels'][group_sizes.index(group_size)]
                    array_storing_info_for_plot.append(data_to_store)
                    x_ticks.append(key)

            plt.subplot(number_of_rows,number_of_columns,i+1)
            plt.boxplot(array_storing_info_for_plot);

            plt.ylabel('Received level of detected echoes, \n dB SPL re 20$\mu Pa$',
                                                    fontsize=10)
            plt.xticks(np.arange(1,len(x_ticks)+1), x_ticks, fontsize=10)
            plt.yticks(np.arange(20,86,6),np.arange(20,86,6), fontsize=10)
            plt.text(1,21,'Hearing threshold', fontsize=10)
            plt.xlabel('d distance (multiples of minspacing)', fontsize=10)
            plt.legend()
            plt.title(f"group size={group_size}")
            print(key_for_dict)
            i+=1
    label= "/received_level_dB/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_radial.png")
    plt.clf() 

    i=0; 
    number_of_columns=3; number_of_rows= int(np.ceil(len(sublist_of_simulation_data_intransect)/number_of_columns))
    plt.figure(figsize= (25,25))
    for simulation_data in sublist_of_simulation_data_intransect:
            groupsizes= group_sizes
            key_for_dict= simulation_data["focal_bat"][0][0]*sign_based_on_distance_from_centermost_point[i]
            nbrs_detected_data= store_probability_of_detection_n_neighbour[key_for_dict]

            plt.subplot(number_of_rows,number_of_columns,i+1)
            plt.boxplot(nbrs_detected_data['neighbours_azimuth'][:], label=f"r={key_for_dict}")
            plt.plot(np.arange(1,7), nbrs_detected_data['median_neighbours_azimuth'], '-*', linewidth=6,
                            alpha=0.5, label='median neighbour detection angle')
            plt.vlines(range(1,len(nbrs_detected_data['5%ile_neighbours_azimuth'])+1),nbrs_detected_data['5%ile_neighbours_azimuth'],
                    nbrs_detected_data['95%ile_neighbours_azimuth'], label='95% data interval')
            #plt.plot(np.arange(1,8), nbrs_detected_data['95%ile_neighbours_azimuth'], 'g')
            #plt.plot(np.arange(1,8), nbrs_detected_data['5%ile_neighbours_azimuth'], 'g', label='95%  neighbour detection range')

            plt.yticks(np.arange(-180,210,30), np.arange(-180,210,30), fontsize=10)
            plt.xticks(np.arange(1,7), group_sizes, fontsize=10)
            plt.legend()
            plt.ylabel('Neighbour detection angle, \n  Azimuth, $^{\circ}$', fontsize=10)
            plt.xlabel('Group size', fontsize=10)


            print(key_for_dict)
            i+=1
    label= "/neighbour_azimuth/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_groupsize.png")
    plt.clf() 


    i=0; 
    d_distances_along_transect= [-6,-5,4,-3,-2,-1,1,2,3,4,5,6]
    number_of_columns=3; number_of_rows= int(np.ceil(len(group_sizes)/number_of_columns))
    plt.figure(figsize= (21,21))
    for group_size in group_sizes:
            array_storing_info_for_plot= []
            iterant_sign=-1
            x_ticks=[]
            for key in d_distances_along_transect:
                    iterant_sign+=1
                    nbrs_detected_data= store_probability_of_detection_n_neighbour[key]
                    data_to_store= nbrs_detected_data['neighbours_azimuth'][group_sizes.index(group_size)]
                    array_storing_info_for_plot.append(data_to_store)
                    x_ticks.append(key)

            plt.subplot(number_of_rows,number_of_columns,i+1)
            plt.boxplot(array_storing_info_for_plot);

            plt.yticks(np.arange(-180,210,30), np.arange(-180,210,30), fontsize=10)
            plt.xticks(np.arange(1,len(x_ticks)+1), x_ticks, fontsize=10)

            plt.ylabel('Neighbour detection angle, \n  Azimuth, $^{\circ}$', fontsize=10)
            plt.xlabel('d distance (multiples of minspacing)', fontsize=10)
            plt.title(f"For group size {group_size}")
            plt.legend()
            print(key_for_dict)
            i+=1
    label= "/neighbour_azimuth/"
    make_dir(dir_to_store_plots+label)
    plt.savefig(dir_to_store_plots+label+"/paired_by_radial.png")
    plt.clf() 