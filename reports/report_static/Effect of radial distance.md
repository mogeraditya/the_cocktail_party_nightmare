#thesis #report #static_simulation

#### 28h july 2025
(side notes : ignore -2 in groupsize 10 some issue is there while loading the data)
#### pipeline of varying radial distance; what is implemented and how?
- to change the position of the focal bat we need the following inputs
	- radial distance ($d$; $d \in \mathbb{N} + \{0\}$): the distance of the new focal bat from the center-most bat. this is entered as integer multiples of the minimum spacing between points ($min\_spacing$). minimum spacing is fixed because all the points are placed by using a Poisson disk sample. 
	- theta ($\theta$; $\theta \in [0,360]$): angle of the focal bat with respect to the x,y coordinates with the center-most bat as the origin. 
- first thing is to ensure a point exists given some radial distance and angle, the following are the criteria we used;
	1. The points need to be in the same quadrant as $\theta$
	2. two lines are drawn based on theta and a preset threshold (0.1 in all our analysis).
		- $y= tan(\theta)*x \pm threshold* min\_spacing$  are our two lines
		- the first criteria is satisfied if a point exists within this bound
	3. if distance from the origin of points selected in 1 and 2 is within the following bound, the point is selected.
		- bound: $min\_spacing*(d \pm threshold)$ 
	4. one point randomly sampled from points that are sublisted from the above selection
	5. if no such point exists, the bats are randomly placed once again and this is carried out until an arrangement with a point satisfying the criteria is found. the search is terminated if no such arrangement is found for 1000 iterations.
- we took transect long the y axis; positive values of $d$ mean $\theta=90\degree$ and negative values of $d$ mean $\theta=270\degree$

Results
1. effect on number of neighbours detected 
	![[Pasted image 20250730163315.png]]
	![[Pasted image 20250730163324.png]]
	
2. effect on probability of detecting atleast n neighbours per call
	![[Pasted image 20250730165140.png]]
	many important points
	- generally, the bats at the back have better probabilty of detection across plots than the centermost bat. the bats at the front have worse probability of detection that the bats in the center. This difference is only exxagerated with increase n (in atleast n neighbours detected) [CHANGE] the exaggeration is visible because of the error bars
	- when the bat is in the frontedge of the group, the difference is the highest, and this difference is less drastic with increase in group size mostly becuase overall predictions drop.
	- the back edge does slightly worse than the bat 1 din front of it consistently. unsure about why, *discuss with TB* for now. this is not there for -2 edge case for group size 10 because of slight error in the data; unsure about where this arises from. 
	- centermost seems to be the worst at detecting atleast one neighbour. [*discuss with TB*] 
	- edge case vs group size (1 for 5; 2 for 10; 3 for 20; 4 for 30; 5 for 50; 6 for 75)
	- front edge always has the worst probability and back edge doesnt always have the best probabilities
3. effect on azimuth of detected neighbors
	![[Pasted image 20250730165316.png]]
	![[Pasted image 20250730172405.png]]
	the bat at the back of the emerging group obviously has more bats in the front that a bat in the front of the group as. therefore as you move from the back edge of the group to the front of the group you expect to go from a uni modal distribution to a more bimodal distribution 
4. effect on received levels of echoes
	![[Pasted image 20250730165347.png]]
	