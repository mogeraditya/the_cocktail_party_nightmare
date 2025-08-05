
#thesis #report #static_simulation 
#### 28th July 2025

There are three main classes of group structure on the basis of heading variation
1. emergence: This situation tries to replicate the acoustic scene when bats are emerging out of a cave. all bats face one direction with a variance of $\pm \theta; \theta \in [0,90]$. This theta was set to $10$ in the static simulations paper. 
2. swarming: this situation tries to replicate the acoustic scene when bats are in a swarming formation. Bats have random heading directions. It is the same as the emergence situation if theta is set to $180$ degrees.
3. milling: this situation tries to replicate the acoustic scene when bats are in a milling formation. This is when bats are flying around in a circular formation. 
### Implementation of the milling formation;
1. pick an arbitrary center as the center around which all the points are spinning around. In our case this is the mean point (mean(all_x_coordinates)), mean(all_y_coordinates))
2. Make multiple sections of the 2D space centered around this point. the number of sections will dictate the smoothness of the concentric circles. For all of the analysis 32 sections are used. that means that the angle between each section is 
    $\theta = \frac{360}{32} = 11.25 \degree$
3. All points within the same sections have the same heading direction. Any given point's heading direction is perpendicular to the border between its section and the very next section (sections are ordered in counter clockwise fashion). For example, for 32 sections, all points within angles 0 to 11.25 degrees are considered as part of one section, from 11.25 to 22.50 degrees are part of the next section. points in the 0 to 11.25 section are facing the direction perpendicular to the border, perpendicular to the line 
	$y= tan(11.25\degree) * x$ 
4. an example of this is given below for multiple number of sections;
	1. ![[Pasted image 20250728162107.png]]
	2. ![[Pasted image 20250728162113.png]]
	3. ![[Pasted image 20250728162118.png]]
	4. ![[Pasted image 20250728162125.png]]

### Parameters;
1. We have used 32 sections for all simulation runs.
2. For the emergence case we looked at 4 different values of $\theta$, [0,10,30,50]
3. We looked at group sizes ranging from [5,10,20,30,50,75]
4. Each combination of these parameters were simulated 100 times. 

### Results
1. effect on number of neighbours detected.
![[Pasted image 20250728170413.png]]
![[Pasted image 20250728170358.png]]
There seems to be no significant differences between different heading variations in the number of neightbours detected.

2. effect on detected neighbour azimuth
![[Pasted image 20250728172929.png]]
3. effect on spl of detected echoes
![[Pasted image 20250728180140.png]]

These dont have significant differences in the response varaibles across the different heading var parameters because the focal bat is the centermost bat.
the auditory scene is not very different between individuals when the focal bat is fixed at the center even when the heading variaiotns are cahnged because on average you have the same set of individuals facing away / towards you. milling case is slightly different, but this difference fails to show up anywhere except in the response in received level of echoes.

But one thing does have difference which is super exciting 
4. effect on probability of detecting atleast n neighbours
	![[Pasted image 20250730130151.png]]
milling sort of forms an upper bound on the >= n neighbours detected per call, i.e., the center bat in the milling case has the highest neighbour detection probabilities per call. This is also reflected in the probability of detecting zero neighbours per call plot, where milling is the lowest among all other formations. 
the one that has the lowest is the emergence case with the highest variation/
Let's subsample to get a distributiion of probabilities.
