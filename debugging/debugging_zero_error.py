import sys
sys.path.append("../simulations/")
from run_simulations import run_multiple_simulations


run_multiple_simulations(name="focal_bat_d", info="trail run",
                                parameter_file_format="./effect_of_position_focal_bat_d/paramsets/*",
                            num_CPUs=4,
                            dest_folder="./effect_of_position_focal_bat_d/store_results_clean/", varying_param= "change_focal_bat_d")   