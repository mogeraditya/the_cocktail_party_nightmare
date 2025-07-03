# The Cocktail Party Nightmare
![](docs_imgs/theCPN_logo.jpg)
#### Author: Thejasvi Beleyur, Max Planck Institute for Ornithology, Seewiesen.
#### Last Updated: November 2019

#### Code associated with the paper:
#### Modelling active sensing reveals echo detection even in large groups of bats, Beleyur & Goerlitz 2019, Proceedings of the National Academy of Sciences

### What is the repository about?
The code in this repository simulate what a bat flying in a group of other bats may be experiencing. Echolocating bats fly in the night and live in caves where they cannot 'see' objects with their eyes like we do. They emit loud calls and listen to the returning echoes to detect their surroundings. See [here](https://www.cell.com/current-biology/fulltext/S0960-9822(05)00686-X) for a quick review. 

Echolocation works very well when individual bats are flying around alone while hunting or commuting. However, when there are  loud sounds, such as those from other bats - the faint returning echoes may not be heard. This means that when a bat flies close to many other bats - it may be flying metaphorically 'blind'. However, many bat species are very social, and live in large groups, fly around together in caves and even emerge in the millions - the Mexican Free-Tailed bat is a classic example. 

The problem of listening to a target signal in the presence of louder sounds that affect signal detection has been called the 'cocktail party problem' by [Cherry 1953](https://asa.scitation.org/doi/10.1121/1.1907229), (also the [wikipedia article](https://en.wikipedia.org/wiki/Cocktail_party_effect)). In the context of bat echolocation, given the orders of magnitude louder bat calls and the large number of calls and echoes that need to be dealt with when flying in a group - [Ulanovsky & Moss 2008](https://www.pnas.org/content/105/25/8491.short) termed group echolocation a 'cocktail party nightmare'. 

This repository simulates the auditory detection of echoes, sound propagation and other relevant phenomena that occur in the cocktail party nightmare and quantifies how many neighbours a bat may be detecting as it flies in a group. For more detail on the simulation implementation please refer to the associated paper <LINK TO PREPRINT HERE>. 


### Requirements:
The code runs on a range of Python 3.10 versions and should run on Linux and Windows systems. All simulation results generated for the paper were done on a Ubuntu 18.04.3 LTS virtual machine. The code was also tested and developed on a Windows 7 system.
To replicate the development environment exactly - utilize the cocktail_3_env.yml file to replicate the exact versions of packages into your environment.

### Installing the required Python environment:
It is a good idea to always use environments for each project. Environments are unique ''boxes' that have their own installation of Python and dependent packages ([article on environments](https://protostar.space/why-you-need-python-environments-and-how-to-manage-them-with-conda)).

I used an [Anaconda](https://docs.anaconda.com/anaconda/) Python installation. (add steps to recreate)


## Running simulations described in the paper: 4 steps to re-running it on your system
All simulation code is in the 'simulations' folder. Every set of simulations is defined by a parameter file which sets the number of simulation runs, call properties, group size and other relevant parameters in the cocktail party nightmare. The simulations themselves are run by calling the ```run_simulations``` module with arguments through a command line call. Let's go through it step by step.

Before you begin trying to run the simulations - first make sure you have  a Terminal window (Unix) or a Command Prompt (Windows) open and move into the 'the_cocktail_party_nightmare' repository.

#### Step 1 : Generate the common parameter file
Each of the simulation scenarios in the 'simulations' folder rely on a common set of parameters that they then alter according to the situation. Generate the common parameter set first to allow the creation of more specific simulation parameter sets. Generate the common simulation parameter file by running the following command 
```
cd simulations/
python create_and_save_common_simparameters.py
```
Running this will lead to the creation of a ```'common_simulation_parameters.paramset'``` file in the folder. The ```.paramset```  extension defines a parameter file. 

##### What is a parameter file?
A parameter file consists of a dictionary with simulation parameters such as group size, call duration, source level, inter-bat spacing etc. 
Having setup the parameters you can simply run the ```make_<simulation_type>.py```' file and to get a bunch of ```.paramset``` files in your folder. The ```make_<simulation_type>.py``` modules in each simulation scenario ('effect_of_group_size, multivariable_simualtions') are the reference examples. 

#### Step 2 : Setting up a simulation by creating specific parameter files
If you want to start by recreating the simulations that have been run in the paper, you can directly run the ```make_<simulation_type>.py``` module directly through an IDE like Spyder or run it from the command line with the following commands. 

Alternatively, if you would like to explore your own parameters, alter parameters in the 'make_params' files through the text editor of your choice. For more documentation on the parameters that can be varied and their data types/formats please see the 'run_CPN' function in the 'the_cocktail_party_nightmare.py' located in the 'CPN' folder. In addition to the documentation, the ```make_<simulation_type>.py``` modules are also useful to understand how the parameters are fed into the simulations. 

At least in the beginning keep the group sizes simulated small (all_group_sizes = [10, 50]), and run only a few simulation runs (number_of_simulation_runs = 1). The code could take a long time to run at larger group sizes depending on your system specs, and it might be nice to get the results quickly in the pilot runs!
```
# move into the directory where the make .py module is. 
cd simulations/effect_of_group_size
python make_params_groupsize.py
```

You should now see multiple ```.paramset``` files in the folder, each of which has the parameters required to initiate simulations for the different group sizes. 

#### Step 3 : Starting the simulations
The simulations are initiated through the ``` run_simulations.py``` module in the simulations folder. The ```run_simulations``` module initiates simulations runs as described by each parameter file in a folder. The simulation outputs are saved into the user-given destination folder. 

#### Starting simulations runs 
```
#Change directories back to the 'simulations' folder. 
cd the_cocktail_party_nightmare/simulations/
# run effect of group size simulations 
python run_simulations -name "groupsize_effect" -info "a trial run to see if things work okay" -param_file "effect_of_group_size/*.paramset" -dest_folder "effect_of_group_size/" -numCPUS 2
```

The ```-name``` argument will set the prefix for all simulation output files. Here all outputs will start with 'groupsize_effect'.

The ```-info``` argument will be saved into the simulation output so you can quickly refer to what it is that was being run without actually loading the information. 

The ```-param_file``` argument refers to the folder and file format of the parameter files to be used to initiate the simulations. 

The ```-dest_folder``` argument refers to the destination folder for all the simulation outputs. 

The ```-numCPUS``` refers to the number of CPUs that are to be used to run the simulations in parallel. Simulations with larger group sizes (>50 bats) can be intensive and noticeable if other programs are running alongside the simulations. It may then be wise to limit the number of CPUs that are being used for the simulations. If the ```-numCPUs``` parameter is not defined the value defaults to the total number available on the device if not specified (typically 4 for your average laptop/PC). 

#### What a succesful simulation initiation looks like : 
![](docs_imgs/succesful_run.png)
You should see the initiation of many parallel simulation runs.

#### Step 4 : The simulation outputs
Every simulation run produces its own ```.simresults``` file. Each simresults file has two objects within it: 

* *simulation_identifiers* : the parameters that were used to run the simulation and the random seed that can be used to re-run the exact same simulation run on any operating system if need be. These include the ```info``` tag, the [uuid](https://en.wikipedia.org/wiki/Universally_unique_identifier), parameter values. 

* *simulation_data* : the results (data) of the simulation run

The *simulation_identifiers* allow each simulation to be uniquely identified by its , and exactly replicated through its random seed.  The  *simulation_data* is the actual data that are used in the analyses.

#### Loading, analysing and handling simulation outputs :
The Jupyter Notebooks in the 'analysis' folder provide a more detailed glimpse of how to proceed with handling and analysing the simulation results - here is a brief glimpse. 

#### Checking out the simulation identifiers
The simulation identifier object is a dictionary with all the relevant identifiers a simulation run has. 
```
# open a Python session or a Jupyter Notebook

# load the .simresult files
import dill 
import numpy as np 
import sys 
sys.path.append('../CPN/') # I'm not sure why but the data loading throws an 
# error if the 'CPN' folder is not in the search path - help appreciated!

# load the simresult file 
with open('<your simulation run output file here>.simresults', 'rb') as simfile:
    one_sim_run = dill.load(simfile)

# unpack the simulation identifiers and simulation data 
sim_id, sim_data = one_sim_run 

# display some of the simulation identifiers
print(sim_id)
```
![](docs_imgs/sim_id_output.png)

#### Checking out the simulation data :
The simulation data is a list with 3 objects in it:  the echoes heard, the sounds in the interpulse interval, group geometry.

```
# unpack the simulation data into its three component objects
echoes_heard, sounds_in_ipi, group_geometry = sim_data 
```

* i) echoes heard : In each simulation the focal bat could hear up to group size -1 echoes bouncing off each of its neighbours. The ```echoes_heard``` object is a numpy array with 0's and 1's. The 1's represent a neighbour detection, the 0's represent a neighbour that was not detected. 

    The ```echoes_heard``` object has group size -1 entries because the focal bat only hears echoes that reflect off *other* bats. The identity of un/detected neighbours is obtained by looking at the index of the ```echoes_heard``` array. 

    If an array were to look like:
    [0,1,0,0,0] - this means the i=1 st index neighbour was detected. To get the xy and heading position of the +1 must be added - and so   the xy position and heading of the i =2 nd bat must be accessed in the group geometry object.

    An example from a simulation run with 50 bats, and thus only 49 entries : 
    ![](docs_imgs/echoes_heard.png)

* ii) sounds in the interpulse interval :  a dictionary with 3 Pandas DataFrames describing the details of each sound's 
path to the focal bat, its time of arrival in the interpulse interval, the distance it travelled, the angle at which it arrived etc. 

  * 1. target echoes :  the echoes generated by the focal bat that are reflected off its neighbours. 
  * 2. secondary echoes : the echoes *not* generated by the focal bat, but from emissions of neighbouring bats.
  * 3. conspecific calls : the neighbouring bat calls. 
  
    Here is an example of the ```target_echoes``` DataFrame detailing the route the echo took as it went from source to received amongst other things. The numbers in the route correspond to the index numbers identifying each bat. If the route is (0,3,0), it corresponds to an echo that was emitted by the focal bat, went to bat 3 and returned to the focal bat. 
    ![](docs_imgs/target_echoes.png)

 *  iii) group geometry : a dictionary with the heading directions and xy positions of all bats in a group. Headings refer to the direction a bat is flying in and aiming its sonar beam at. 
 
       orientations : This object has the direction in which each bat was aiming its beam. The angles are between 0-360 degrees, with the 0 degrees being at 3 o'clock and increasing in counter-clockwise fashion. The orientations object is a Nbats x 1 np.array. The 0th index orientation corresponds to the focal bat. 
          
    ![](docs_imgs/orientations.png)
    
       positions : The xy co-ordinates of all bats. The 0th row corresponds to the xy positions of the focal bat.
                      Nbats x 2 np.array with xy coordinates of the group of bats.
    ![](docs_imgs/positions.png)

#### Performing further analyses of the already run results:
The simulation results generated add up to a few hundred GB and a lot of individual files and so I have not yet been able to compress and upload them onto a public archive. If you are interested in analysing the raw data I would be happy to share the raw data through a suitable medium/channel! 

### Associated experimental data and code-base
This repository by itself only has the code and basic associated simulation parameter data to get the simulations running. The simulation parameters are based on previously published literature and two experiments specifically performed to parametrise the simulations. The bistatic and monostatic target strengths of bats were measured along with how bats 'shadow' or block sound. Head to this [Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3469845.svg) repository to read about the experimental methodology, code, raw data and results involved in these experiments. 
If you find the code/data or playback files useful in your work please consider citing the paper associated with the work.

### License:
All code in the_cocktail_party_nightmare is release under an MIT License. See [LICENSE](./LICENSE) for more detail. 

### Citing the cocktail party nightmare code-base
If you find any of the code in this repository useful in your research please cite the [associated paper](https://www.pnas.org/content/early/2019/12/09/1821722116):
Beleyur T., Goerlitz, H.R. Modelling active sensing reveals echo detection even in large groups of bats

In addition, this is the link to the publicly available [preprint](https://www.biorxiv.org/content/10.1101/817734v1):

```Beleyur, T., & Goerlitz, H. R. (2019). Active sensing in groups:(what) do bats hear in the sonar cocktail party nightmare?. BioRxiv, 817734 ```

#### Read till here and want to know more? 
 Please do not hesitate to contact us (thejasvib@gmail.com, hgoerlitz@orn.mpg.de) or raise an Issue on Github if you
 
 a) have tried to replicate the simulations and are having problems
 
 b) found bugs in code
 
 c) have ideas on how this project can be developed further
 
 
