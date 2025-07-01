import collections
import unittest
import numpy as np
import sys
sys.path.append("../CPN/")
from the_cocktail_party_nightmare import generate_line_thresholds_given_theta_d, check_if_given_point_between_two_lines 
from the_cocktail_party_nightmare import check_quadrant, check_distance, generate_subset_of_points_given_d, generate_new_focal_bat_given_angle

class TestingNoncentralAllocation(unittest.TestCase):
    def setUp(self):
        
        self.points_x= [1,-1,1,-1,0]
        self.points_y= [1,1,-1,-1,0]
        self.lines= (0,1)
        self.points_to_test= [[1,0], [1,1], [-1,-1], [-1, 1.2], [1, -1.2]]
        self.kwargs= {}
        self.kwargs["Nbats"]= 100
        self.kwargs["min_spacing"]= 0.5

    def test_check_quadrant(self):
        theta_true= 0; theta_false= 3*np.pi/4
        point= np.array(self.points_to_test[0])
        bool_theta_true= check_quadrant(theta_true,point)
        bool_theta_false= check_quadrant(theta_false,point)
        self.assertTrue(bool_theta_true and not bool_theta_false)

    def test_check_if_given_point_between_two_lines(self):

        expected_true_false= [True, True, True, False, False]
        computed_true_false= [check_if_given_point_between_two_lines(self.lines, i) 
                              for i in self.points_to_test]
        # print(computed_true_false)
        self.assertTrue((expected_true_false==computed_true_false))

    def test_check_distance(self):
        number_of_d = 2; theta=np.pi/4
        expected_true_false= [True, False, False, False, False]
        computed_true_false= [check_distance(number_of_d,self.kwargs["min_spacing"], i, theta, self.lines)
                                for i in self.points_to_test]
        print(computed_true_false)
        self.assertTrue((expected_true_false==computed_true_false))
        
    #the other function just implements already tested functions and therefore doesnt require unit test
if __name__ == '__main__':

    unittest.main()