# -*- coding: utf-8 -*-
""" Module to implement acoustic shadowing along with spherical spreading
of sound

Created on Mon Jun 17 16:08:50 2019

@author: tbeleyur
"""
import time
import sys 
sys.path.append('..//bridson//')
sys.path.append('..//')
import numpy as np 
import pandas as pd
import scipy.spatial as spatial
import statsmodels.api as sm


def soundprop_w_acoustic_shadowing(start_point, end_point, all_other_points,
                                   **kwargs):
    '''Calculate the received level of a sound emitter at the start point 
    and reaching the end point after potentially passing through other points on 
    the way. 
    
    Each point sound passes through creates a drop in the intensity because of 
    acoustic shadowing. 

        

    Parameters
    ----------

    start_point : 1 x 2 array like

    end_point :  1 x 2 array like

    all_other_points : Nbats-2 x 2 array like
                       xy coordinates of all points between start and end point
                       in a rectangular area of given width 
      
    Keyword Arguments
    ----------------
    implement_shadowing : Boolean. If True then shadowing calculations are done, else only simple spherical spreading. 

    R : float. straight line distance between soure and receiver

    rectangle_width : float >0. width of the rectangle 

    emitted_source_level : dictionary with key 'dBSPL' and 'ref_distance' ,
                    indicating source level in dB SPL re 20muPa and reference distance in metres. 
    
    R : float >0. straight line distance between source and receiver. This is used
        in case there are no obstacles between them. 
   
    acoustic_shadowing_model : statsmodel object that allows calculation of how much
                               shadowing will be observed.      
                               
    min_spacing : float> 0. 


    Returns
    -------
    received_level : float. 
                    received level of the sound 
    '''    
    
    if kwargs.get('implement_shadowing'):
        all_points_between = get_points_in_between(start_point, end_point, 
                                                       all_other_points, **kwargs)
        
        received_level = calc_RL(kwargs['R'],kwargs['emitted_source_level']['dBSPL'],
                                             kwargs['emitted_source_level']['ref_distance'])
        
        num_obstacles = all_points_between.shape[0]        
        if num_obstacles >= 1:
            acoustic_shadowing = calculate_acoustic_shadowing(num_obstacles, 
                                                              **kwargs)
            received_level += acoustic_shadowing          
    else:
        received_level = calc_RL(kwargs['R'],kwargs['emitted_source_level']['dBSPL'],
                                             kwargs['emitted_source_level']['ref_distance'])
    return(received_level) 



def get_distances_between_points(xy_between, start, end):
    '''
    
    '''
    all_xy =  np.row_stack((start, xy_between, end))
    distance_matrix = spatial.distance_matrix(all_xy, all_xy)
    distance_to_source = np.argsort(distance_matrix[:,0])

    points_sorted = all_xy[distance_to_source,:]
    distances_sorted = spatial.distance_matrix(points_sorted, points_sorted)

    num_distances = all_xy.shape[0] - 1
    point_2_point_distances = np.zeros(num_distances)
    for i, (point0, point1) in enumerate(zip(range(all_xy.shape[0]-1),
                                         range(1,all_xy.shape[0]))):
        point_2_point_distances[i] = distances_sorted[point0, point1]
    return(point_2_point_distances)
    

def calculate_acoustic_shadowing(num_obstacles,
                                         **kwargs):
    '''Calculates received level of a call with acoustic shadowing included. 
    The received level of the call with shadowing is calculated with an iterative 
    application of the bistatic sonar equation. The TS used here is the bistatic 
    target strength at emitter-receiver angular separations of 180 degrees. 
    
    Parameters
    ----------
    num_obstacles : int >1 .
                    Number of obstacles between receiver and emitter.


    Keyword Arguments
    -----------------

    acoustic_shadowing_model : statsmodel object
                                A statistical model that allows calculation of
                                the amount of acoustic shadowing in dB. 
                                For predictions the model accepts a 
                                pd.DataFrame which has the following 
                                columns (this might depend on the exact model 
                                loaded too!)
                                
                                obstacles 
                                spacing

    min_spacing : float>0. 
                    Separation between bats/obstacles
                    see the_cocktail_party_nightmare

    Returns
    -------
    shadowing_reduction : float. 
                        Reduction of received level due to shadowing in dB. 
  
    '''
    no_obstacle = pd.DataFrame(data={'obstacles':[0],
                                     'spacing':[kwargs['min_spacing']],
                                     })
    with_obstacles = pd.DataFrame(data={'obstacles':[num_obstacles],
                                     'spacing':[kwargs['min_spacing']],
                                     })
    #convert_to_categorical(no_obstacle, 'spacing')
    #convert_to_categorical(with_obstacles, 'spacing')
    level_w_obstacles = kwargs['acoustic_shadowing_model'].predict(with_obstacles) 
    level_wo_obstacles = kwargs['acoustic_shadowing_model'].predict(no_obstacle) 
    shadowing_reduction = float(level_w_obstacles - level_wo_obstacles)

    return(shadowing_reduction)

def convert_to_categorical(df, column):
    '''
    '''
    df[column] = pd.Categorical(df[column])
    return(df)


def calc_RL(distance, SL, ref_dist, **kwargs):
    '''calculates received level only because of spherical spreading.

    Parameters
    -----------

    distance : float>0. receiver distance from source in metres.

    SL : float. source level in dB SPL re 20 muPa at the reference distance.

    ref_dist : float >0. distance at which source level was measured in metres.
                Typically 1metre by convention.

    Keyword Arguments
    -----------------
    atmospheric_attenuation : float <= 0. 
                              Atmospheric attenuation in dB/m. 
                              This has to be negative number. 
                              Defaults to no atmospheric attenuations (0 dB/m )
    


      Returns
    -------

    RL : received level in dB SPL re 20muPa.

    '''
    RL = SL - 20*np.log10(float(distance/ref_dist))
    RL += kwargs.get('atmospheric_attenuation', 0)*distance

    return(RL)

def get_points_in_between(start_point, end_point, all_other_points,
                          **kwargs):
    '''
    
    Parameters
    ----------
    start_point : 1x2 array like
                  xy coordinates of starting point  

    end_point : 1x2 array like
                xy coordinates of end points

    all_other_points : N x 2 array like
                       xy coordinates of all other points  
    
    Keyword Arguments
    -----------------
    rectangle_width : float >0.
                      The width of the rectangle between the start and end point.

    Returns
    -------
    points_between : Mpoints_between x 2 np.array where Mpoints can be >= 0. 

    '''
    rectangle_limits, rotation_matrix = make_rectangle_between_2_points(start_point, 
                                                       end_point,**kwargs)

    points_between = get_points_in_rectangle(rectangle_limits, start_point,
                                                   all_other_points, rotation_matrix)

    return(points_between)

def get_points_in_between_thecircleversion(start_point, end_point, 
                                          all_other_points,**kwargs):
    '''Take 2 at getting perhaps a faster version of the 
    previous get_points_in_between function. 
    
    It is fast *and* dirty ... and doesn't quite apply when many bats are packed tightly
    together ... as long as the 'rectangle_width' is decently large -- then it
    should be okay..

    
    '''
    # get line equation from A to B
    diff_x_y = end_point-start_point
    vertical, m = calculate_slope(diff_x_y)
            
    numpoints = 100 # choose a default density for now 
    
    points_along_line = np.zeros((numpoints, 2))
    
    if not vertical:
        points_along_line[:,0] = np.linspace(start_point[0],end_point[0],
                                                 numpoints) # x coordinates
        c = solve_for_intercept(start_point,end_point,m)
        points_along_line[:,1] = m*points_along_line[:,0] + c # y coordinates 
    
    else:
        points_along_line[:,0] = start_point[0] # x coordinates
        points_along_line[:,1] =  np.linspace(start_point[1], end_point[1],
                                 numpoints)# y coordinates 
    
    # get the distance from each of the points to all other points
    distance_from_line = spatial.distance_matrix(points_along_line, 
                                                         all_other_points)
    within_r_dist_from_line = distance_from_line<= kwargs['rectangle_width']

    point_ids = np.argwhere(np.any(within_r_dist_from_line, axis=0))
    return(all_other_points[point_ids])
        
        

def calculate_slope(deltax_deltay):
    '''
    Parameters
    ----------
    deltax_deltay : 1 x array-like. [X2-X2, Y2-Y1 ]

    Returns
    -------
    vertical : boolean. True if the line is vertically oriented (deltaX==0)

    slope : float. 
    
    '''
    zeros_present = tuple((deltax_deltay == 0).tolist())

    if sum(zeros_present)>0:
        return(slope_dict[zeros_present])
    
    else:
        return(False, deltax_deltay[1]/float(deltax_deltay[0]))

slope_dict = {(True,False) : (True, np.nan) ,# xchanges y doesnt
              (False, True): (False, 0.0), #y doesnt, x changes
              (True, True): Exception('Zero slopes for both not possible!!')}
        
     
    

def solve_for_intercept(x1y1, x2y2,m):
    '''
    '''
    x1,y1 = x1y1
    x2,y2 = x2y2
    
    c = (y1 + y2 - m*x1 - m*x2)/2.0
    return(c)
    
    

    

def make_rectangle_between_2_points(A, B, **kwargs):
    '''First calculate the relative coordinates of B wref to A. 
    Then draw a straight line between 0,0 and B_rel and 'undo' the
    slope. To this vertical line not apply rectangular bands on left and right.
    
    Output bottom left and top right vertices along with the rotation matrix
    along with the rotation matrix used to undo the slope of the 0,0-B_rel line
    for application to other poitns. 

    Parameters
    ----------
    A, B : 1x2 array like. xy coordinates of start(A) and end(B) points

    Keyword Arguments
    -----------------
    rectangle_width : float>0. The width of the rectangle between A and B. 

    Returns
    -------
    corner_limits : tuple.
                    Consisting of x0,x1,y0,y1

    '''
    # treat A as origin, calculate slope between B and A
    B_rel = B-A 
    # 'un-rotate' B and thus form a vertical rectangle easily
    theta = np.arctan2(B_rel[1], B_rel[0])

    #theta_tobe_rotated = np.remainder(theta, np.pi/2)
    theta_tobe_rotated = np.pi/2.0 - theta
    rotation_matrix = rot_mat(theta_tobe_rotated)
    B_rotated = np.dot(rotation_matrix, B_rel)
    x0, x1 = -kwargs['rectangle_width']*0.5, kwargs['rectangle_width']*0.5
    y0, y1 = 0, B_rotated[1]

    return([x0,x1,y0,y1], rotation_matrix)
	
def get_points_in_rectangle(corner_limits, startpt,
                            many_points, rotn_matrix):
    ''' Many points are checked if they are within 
    the rectangle defined by the bottom left point (x0,y0) and 
    the top right corner (x1,y1).

    Corner limits is a tuple with four entries (x0,x1,y0,y1)
    x0,x1 the x coordinates defining the width of the rectangle
    y0,y1 the y coordinates defininf the height of the rectanlge
   
    rotn_matrix is the rotation matrix to rotate the many points into 
    the same frame of reference as the rectangle. it is in the 
    form of a 2 x 2 array with the form described 
    [here](https://en.wikipedia.org/wiki/Rotation_matrix)



    '''
    x0,x1,y0,y1 = corner_limits
    relative_posns = many_points - startpt
    rotated_pts = np.apply_along_axis(dot_product_for_rows, 1, relative_posns,
                                      rotn_matrix)
    
    within_x = np.logical_and(rotated_pts[:,0] >= np.min([x0,x1]),
                              rotated_pts[:,0] <= np.max([x1,x0]))

    within_y = np.logical_and(rotated_pts[:,1] >= np.min([y0,y1]),
                              rotated_pts[:,1] <= np.max([y0,y1]))
    
    within_pts = np.logical_and(within_x, within_y)
    return(many_points[within_pts])

def dot_product_for_rows(xy_row, rotation_matrix):
    return(np.dot(rotation_matrix, xy_row))

def dot_product_w_sum(xy_row, rotation_matrix):
    return(np.sum(rotation_matrix*xy_row, 1))


def rot_mat(theta):
    rotation_matrix = np.float32(np.row_stack(([np.cos(theta), -np.sin(theta)],
                                    [np.sin(theta), np.cos(theta)])))
    return(rotation_matrix)

# import pickle
# x= r'/home/adityamoger/Documents/GitHub/cocktail_clone/data/acoustic_shadowing_model.pkl'
# with open(x, 'rb') as file:
#     data = pickle.load(file)

# with open("df3.pkl", "wb") as f:
#     pickle.dump(data, f, encoding="bytes")

# print(data)
  
if __name__ == '__main__':
#    kwargs = {'rectangle_width':0.1, 'implement_shadowing':True,
#              'emitted_source_level': {'dBSPL':90, 'ref_distance':1.0}}
#    kwargs['shadow_TS'] = [-15]
#    np.random.seed(82319)
#    #otherpts = np.random.normal(0,5,2000).reshape(-1,2)
#    
#    #otherpts = np.array(([1,0],[1,0.05]))
#    #print(get_points_in_between(np.array([2,0]), np.array([0,0]), otherpts, **kwargs ) )
#    start = time.time()
#    numpoints = [5,90, 100, 1000, 2000, 4000, 8000, 10000, 100000]
#    for num_points in numpoints:
#        y_coods = np.random.choice(np.arange(0.1, 5, 0.01),num_points)
#        x_coods = np.tile(0.05,num_points)
#            
#        between_points = np.column_stack((x_coods, y_coods))
#        q=get_points_in_between(np.array([0,0]), np.array([0,10]), between_points,
#                                   **kwargs)
#    print(time.time()-start)
    kwargs = {}
    kwargs['bats_xy'] = np.array(([0,0],
                                  [1,0],
                                  [2,0]))
    kwargs['focal_bat'] = np.array([0,0])
    kwargs['R'] = 2.0
    kwargs['implement_shadowing'] = False
    kwargs['rectangle_width'] = 0.3
    kwargs['acoustic_shadowing_model'] = sm.load('/home/adityamoger/Documents/GitHub/cocktail_clone/data/acoustic_shadowing_model_p3.pkl')
    kwargs['min_spacing'] = 1.0
    kwargs['emitted_source_level'] = {'dBSPL':90, 'ref_distance':1.0}

    A = soundprop_w_acoustic_shadowing( kwargs['focal_bat'],
                                    kwargs['bats_xy'][-1,:],
                                    kwargs['bats_xy'],
                                    **kwargs)
    print(A)