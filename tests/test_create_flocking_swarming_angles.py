import os 
import glob
import sys
sys.path.append("../CPN/")
from create_milling_swarming_angles import generate_heading_vectors_concentric_circles, generate_heading_vectors_random_all_angles

import unittest
import numpy as np
# from scipy.stats import chisqaured

class Testing_creating_milling_swarming_angles(unittest.TestCase):
    def setUp(self):
        
        self.points_x= [1,-1,1,-1,0]
        self.points_y= [1,1,-1,-1,0]
        self.kwargs= {}
        self.kwargs["Nbats"]= 100

    def test_concentric_circle_heading_vectors(self):
       
        nearby_points= np.column_stack((self.points_x, self.points_y))
        expected_heading_vector= 180
        computed_heading_vector= generate_heading_vectors_concentric_circles(nearby_points, 4)[0]
        vector_match = np.array_equal(computed_heading_vector, expected_heading_vector)
        self.assertTrue(vector_match)

    # def test_random_vectors_all_angles(self):
    #     nearby_points= np.column_stack((self.points_x, self.points_y))
    #     computed_heading_vectors= generate_heading_vectors_random_all_angles(self.kwargs)
    #     chi_squared_test= 
    #     self.assertTrue(vector_match)


if __name__ == '__main__':

    unittest.main()