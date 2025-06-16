import collections
import unittest
import numpy as np
import sys
sys.path.append("../CPN/")
from the_cocktail_party_nightmare import choose_a_noncentral_bat

class TestingNoncentralAllocation(unittest.TestCase):
    def setUp(self):
        self.points_x= [1,-1,1,-1,0]
        self.points_y= [1,1,-1,-1,0]

    def test_focal_bat_change_theta(self):
        nearby_points= np.column_stack((self.points_x, self.points_y))
        centermost_point= [0,0]
        r,theta= np.sqrt(2), np.pi/3; non_central_bat=(r,theta)
        expected_focal_bat= np.array([1,1])
        x, calculated_focal_bat= choose_a_noncentral_bat(non_central_bat, nearby_points, centermost_point)
        centremostpt_match = np.array_equal(calculated_focal_bat, expected_focal_bat)
        self.assertTrue(centremostpt_match)

    def test_focal_bat_change_r(self):
        nearby_points= np.column_stack((self.points_x, self.points_y))
        centermost_point= [0,0]
        r,theta= 2, np.pi/4; non_central_bat=(r,theta)
        expected_focal_bat= np.array([1,1])
        x, calculated_focal_bat= choose_a_noncentral_bat(non_central_bat, nearby_points, centermost_point)
        centremostpt_match = np.array_equal(calculated_focal_bat, expected_focal_bat)
        self.assertTrue(centremostpt_match)


if __name__ == '__main__':

    unittest.main()