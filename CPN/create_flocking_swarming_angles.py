import numpy as np
import math
import matplotlib.pyplot as plt

def flocking(set_of_all_points, resolution):
    '''Goal is to create heading directions in a flock like pattern. 
        * Take the center point and make resolution many equal angled sections of 2*pi. 
        * Assign heading angle based on which section they belong to. 
        * Every points within the same sector have the same heading direction and is perpendicular to the bondary of the next sector. 
        * This function is implemented in an anticlockwise fashion.

     Parameters
    ----------

        set_of_all_points : Npoints x 2 np.array. With X and Y coordinates of points
        resolution: int

      Returns
    -------

        vector_angles : Npoints x 1 np.array.

    '''
    mean_x= np.mean(set_of_all_points[:,0])
    mean_y= np.mean(set_of_all_points[:,1])
    mean= np.array([mean_x, mean_y])
    #this is our center point; now we setup a angle based on resolution
    theta= (2 * np.pi)/ resolution
    transform_points_to_new_center= [i-mean for i in set_of_all_points]
    print(transform_points_to_new_center)
    angles_of_new_points= [np.arctan(i[1]/i[0]) for i in transform_points_to_new_center]

    store_final_vector_angle_id=[]; i=0
    
    for angle in angles_of_new_points:
        point= transform_points_to_new_center[i]
        if point[0]>=0:
            if np.abs(angle)== angle:
                id= np.ceil(angle/theta)
            else:
                id= np.ceil(angle/theta)
        if point[0]<0:
            if np.abs(angle)== angle:
                new_angle=angle-np.pi
                id= np.ceil(new_angle/theta)
            else:
                print(point)
                new_angle=angle+np.pi
                id= np.ceil(new_angle/theta)
        store_final_vector_angle_id.append(id); i+=1
    vector_angles= []
    for id in store_final_vector_angle_id:
        # angle_at_0_id= (np.pi/2)
        vector_angles.append(np.degrees((np.pi/2)+id*theta))
        # else:
    
    print(angles_of_new_points, store_final_vector_angle_id)

    return vector_angles

def swarming(**kwargs):
    '''Goal is to create heading directions in a swarm like pattern. 
        * random heading directions are chosen from -Pi to +Pi.

     Keyword Arguments
    ----------
        Nbats

      Returns
    -------

        headings : Npoints x 1 np.array.

    '''
    min_heading, max_heading = -180, +180
    headings = np.random.choice(np.arange(min_heading, max_heading+1),
                            kwargs['Nbats'])
    return headings

# print(np.arctan(-1))
# print(np.pi/4)
# print(np.ceil((np.pi/6)/(np.pi/4)))
