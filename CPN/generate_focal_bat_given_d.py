import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from the_cocktail_party_nightmare import generate_surroundpoints_w_poissondisksampling
    
#use the poisson disk sample code
#i need to sent the threshold thingy and also implement the inter bat spacing thingy

class No_situation_for_given_nd(Exception):
    """Still an exception raised when uncommon things happen"""
    def __init__(self, message, payload=None):
        self.message = message
        self.payload = payload # you could add more args
    def __str__(self):
        return str(self.message)


def generate_line_thresholds_given_theta_d(theta, d, threshold):
    
    #threshold is in terms fraction of d/ how many times of d is the threshodl
    bound= d*threshold

    #the middle line is at origin; slope is in terms of theta from 0 (radians)
    intercept= bound/np.cos(theta)
    slope= np.tan(theta)

    return (slope, intercept)

def check_if_given_point_between_two_lines(lines, point):
    #lines are parallel
    slope, intercept= lines
    sign_for_upper_line= np.sign(point[1]-slope*point[0]-intercept)
    sign_for_lower_line= np.sign(point[1]-slope*point[0]+intercept)
    if sign_for_upper_line==sign_for_lower_line:
        return False
    # if sign_for_upper_line==0 or sign_for_lower_line==0:
    #     return True
    else: 
        return True
    
def check_quadrant(theta, point):
    #check if point lies in the correct quadrant
    if (point[0], point[1]) == (0,0):
        return True
    if point[0]==0:
        #add a small float
        point[0]+= np.cos(np.pi/2)
    if point[1]==0:
        point[1]+= np.cos(np.pi/2)

    sign_of_sin_cos= (np.sign(np.cos(theta)), np.sign(np.sin(theta)))
    sign_of_point= (np.sign(point[0]), np.sign(point[1]))

    if sign_of_point==sign_of_sin_cos:
        return True
    else:
        return False

def dist_from_origin(point):
    return np.sqrt(point[0]**2 + point[1]**2)

def check_distance(number_of_d,d, point, theta, line_slope_intercept, threshold):
    if check_if_given_point_between_two_lines(line_slope_intercept, point) and check_quadrant(theta, point):
        distance_upper_limit= d*(number_of_d+threshold); distance_lower_limit= d*(number_of_d-threshold)
        distance_from_origin= dist_from_origin(point)
        if distance_from_origin>distance_lower_limit and distance_from_origin<=distance_upper_limit:
            return True
        else:
            return False
    else:
        return False
    
def generate_subset_of_points_given_d(nearby_points, centremost_point, number_of_d, d, theta, threshold):
    transform_nearby_points_st_center_is_origin= nearby_points-centremost_point
    centre_and_other_pts_transformed = np.row_stack((np.array([0,0]), transform_nearby_points_st_center_is_origin))
    # d= kwargs['min_spacing']
    line_slope_intercept= generate_line_thresholds_given_theta_d(theta, d, threshold)    
    subset_of_points= []
    indices= []
    iterant=0
    for point in centre_and_other_pts_transformed[1:]:
        if check_distance(number_of_d, d, point, theta, line_slope_intercept, threshold):
            subset_of_points.append(point)
            indices.append(iterant)
        iterant+=1
    return np.array(subset_of_points), np.array(indices)

def generate_new_focal_bat_given_angle(number_of_d, theta, threshold, **kwargs):

    d= kwargs['min_spacing']
    subset_of_points=[]
    iterant=0
    while len(subset_of_points)==0 and iterant<1000:
        nearby_points,centremost_point = generate_surroundpoints_w_poissondisksampling(kwargs['Nbats'],
                                                                    kwargs['min_spacing'],
                                                                    **kwargs)
        subset_of_points, indices_of_points= generate_subset_of_points_given_d(nearby_points, centremost_point, number_of_d, d, theta, threshold)
        iterant+=1
    
    if iterant==1000:
        raise No_situation_for_given_nd(f"Failure to generate scenario in which there exists no point that is {number_of_d} number of step(s) away from the centermost bat for groupsize of  {group_size}")

    random_index= np.random.randint(0, len(subset_of_points))
    focal_bat_new= subset_of_points[random_index]+ centremost_point; focal_bat_new_index= indices_of_points[random_index]
    nearby_points_without_new_focal_bat= np.delete(nearby_points, focal_bat_new_index, axis=0)
    # print(nearby_points_without_new_focal_bat)
    new_nearby_points= np.row_stack((centremost_point, nearby_points_without_new_focal_bat))
    
    return new_nearby_points, focal_bat_new, nearby_points, centremost_point


# line_slope_intercept= (1, 0.4); point= [1,0.5]; angle=np.pi/3
# if check_if_given_point_between_two_lines(line_slope_intercept, point) and check_quadrant(angle, point):
#     print("slay")
# else:
#     print("cunttty")

group_size = 75
radial_dist= 4
interpulse_duration =  0.1
call_duration =  0.0025
shadowing = True
source_level =  100
spacing = 0.5
group_heading_variation = 10
atmospheric_absorption =  -1
number_of_simulation_runs =  100

simulation_parameters= {}
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
# plt.figure(figsize=(10,25))
# i=0
# while i<10:
#     try:
#         new_nearby, new_focal, nearby, focal= generate_new_focal_bat_given_angle(1, 3*np.pi/2, 0.25, **simulation_parameters)
#     except Exception:
#         print("aiyo amma tai")

#     plt.subplot(5,2,i+1)
#     plt.scatter(new_nearby[:,0], new_nearby[:,1])
#     plt.scatter(new_focal[0], new_focal[1], color="black", alpha= 0.2, label="new focal bats")
#     plt.subplot(5,2,i+2)
#     plt.scatter(nearby[:,0], nearby[:,1])
#     plt.scatter(focal[0], focal[1], color="red", alpha= 0.2, label="centermost bats")
#     i+=2

# plt.show()

# i=0
# while i<10:
#     new_nearby, new_focal, nearby, focal= generate_new_focal_bat_given_angle(2, 3*np.pi/2, 0.1, **simulation_parameters)
#     x= [focal[0], new_focal[0]]
#     y= [focal[1], new_focal[1]]
#     plt.plot(x,y)
#     plt.scatter(focal[0], focal[1])
#     plt.scatter(new_nearby[:,0], new_nearby[:,1], alpha=0.2)
#     i+=1
# plt.show()
# print(new_nearby, new_focal)
import heapq
import scipy.spatial as spl
#find nearest neighbours
def find_n_nearest_points(point, nearby, n):
    centre_and_other_pts = np.row_stack((point, nearby))
    distances_from_centre = spl.distance_matrix(centre_and_other_pts,
                                                centre_and_other_pts)[1:,0]
    smallest= heapq.nsmallest(n, distances_from_centre)
    print(smallest)
    index= [np.where(distances_from_centre==i)[0] for i in smallest]
    points= [nearby[i][0] for i in index]
    print(points)
    return np.array(points)

# new_nearby, new_focal, nearby, focal= generate_new_focal_bat_given_angle(radial_dist, 3*np.pi/2, 0.1, **simulation_parameters)
# nearest_5= find_n_nearest_points(focal, nearby, 3)
# plt.scatter(nearest_5[:,0], nearest_5[:,1])
# plt.scatter(focal[0], focal[1])
# plt.scatter(nearby[:,0], nearby[:,1], alpha=0.2)
# plt.show()
# nearest_5= find_n_nearest_points(new_focal, new_nearby,3)
# plt.scatter(nearest_5[:,0], nearest_5[:,1])
# plt.scatter(new_focal[0], new_focal[1])
# plt.scatter(new_nearby[:,0], new_nearby[:,1], alpha=0.2)
# plt.show()

# i need to look at if all the near neighbours lie within a give angle band lowkey. 
#just draw everything once?
i=0
while i<10:
    new_nearby, new_focal, nearby, focal= generate_new_focal_bat_given_angle(2, 3*np.pi/2, 0.1, **simulation_parameters)
    # nearest_5= find_n_nearest_points(focal, nearby, 3)
    nearest_5_new= find_n_nearest_points(new_focal, new_nearby,3)
    c=np.random.rand(3,)
    x= [focal[0], new_focal[0]]
    y= [focal[1], new_focal[1]]
    for point in nearest_5_new:
        x= [point[0], new_focal[0]]
        y= [point[1], new_focal[1]]
        plt.plot(x,y, c=c)
    plt.scatter(focal[0], focal[1])
    plt.scatter(new_nearby[:,0], new_nearby[:,1], alpha=0.2)
    i+=1
plt.show()