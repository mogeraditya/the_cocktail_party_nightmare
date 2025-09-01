import the_cocktail_party_nightmare as CPN
import matplotlib.pyplot as plt
import numpy as np
from create_milling_swarming_angles import generate_heading_vectors_concentric_circles, generate_heading_vectors_random_all_angles
from generate_focal_bat_given_d import generate_subset_of_points_given_d
import math
np.random.seed(78464)
sim_data={}; sim_data["noncentral_bat"]=False; sim_data["Nbats"]=100
nearby, focal= CPN.generate_surroundpoints_w_poissondisksampling(npoints=100, nbr_distance=0.5, **sim_data)
set_of_all_points= np.row_stack((focal, nearby))
# vector_angles= generate_heading_vectors_concentric_circles(set_of_all_points=set_of_all_points, resolution=32)
# vector_angles = generate_heading_vectors_random_all_angles(**sim_data)
min_heading, max_heading = 0 - 10, 0 + 10
vector_angles = np.random.choice(np.arange(min_heading, max_heading+1),100)
        
points_edge_positive= generate_subset_of_points_given_d(nearby, focal, 2, 0.5, np.pi/2, 0.1)
print(points_edge_positive[0])

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

u_v_s= compute_u_v_for_plotting(set_of_all_points, vector_angles)

plt.figure(figsize= (9,9))
plt.scatter(nearby[:,0], nearby[:,1], color = "b")
plt.scatter(focal[0], focal[1], label="centermost bat", color= "r")
for u_v in u_v_s:
    plt.quiver(*u_v, width = 0.005)
plt.xticks([])
plt.yticks([])
plt.savefig("analysis/plots_for_poster/emergence10.png")
plt.show()

# i=0
# while i<100:
#     nearby, focal= CPN.generate_surroundpoints_w_poissondisksampling(npoints=100, nbr_distance=0.5, **sim_data)
#     set_of_all_points= np.row_stack((focal, nearby))
#     vector_angles= generate_heading_vectors_concentric_circles(set_of_all_points=set_of_all_points, resolution=3)
#     u_v_s= compute_u_v_for_plotting(set_of_all_points, vector_angles)
#     # plt.scatter(focal[0], focal[1])
#     plt.quiver(0,0,*u_v_s[0][2:])
#     i+=1
# plt.show()