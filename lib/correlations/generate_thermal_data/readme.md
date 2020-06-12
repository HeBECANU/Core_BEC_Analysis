# thermal correlation sampling
Bryce Henson

Want to sample (m) counts from a thermal distribution up to some maximum order (p)
probability distrbution(PD)
contitional probability distrbution(CPD)


## rough algo sketch
- starting with a gaussian distrbution sample a single count,using rejection sampling with a proposal density that encloses this distribution
- this then changes the probability distribution by adding a gaussian of the g2 corr. length, and with the amplitude of the thermal distribution at that location
- now readjust the proposal density so that it encloses the probability distrbution
  - will require finding the sum of the thermal PD and all the correlation orders
  - could sample the entire distribution at less than a correlation length( maybe divided by max corr. order)
  - could frame as a search problem starting from each count
  - maybe there is a nice property of gaussians that we can exploit
- now repeat
  - sample a count using rejection sampling 
    - give a sample of x
	- generate a random number between 0 and the proposal density amplitude
	- calcluate the contitional proabablity
	  - thermal 
	  - correlation enghancements
	    - g2 from differences to all other counts
		- g3 ect... up to order p
	- if conditional probablity lies abvove the random number then sucess store that count and the thermap PD amp
	- otherwise repeat
  - readjust the proposal density so that it encloses the conditional proabablity
  
  
  
## some thoughts on complexity
- restrict the interval

for only g2, n counts
- Sampling
  - sample from the prop.dist (c1)
  - calculate the thermal dist (c2)
  - calculate distance differences (n) and coresponding gaussian amp to give conditional proabablity
  - sampling eff int.PD/int. prop.den
- sample the CPD to adjust the prop.dist. n(c2+n)
- repeat unitl have m counts

sum_n=1^m [n+ (n(c2+n))]
to leading term
m^3


for max order p
sum_n=1^m [sum_i=1^p (n^i)  + (n(c2+sum_i=1^p (n^i)))]
to leading order
sum_n=1^m [n^p  + (n(c2+n^p))]
m^(p+1)


for g3 1e3 counts
1e12 operations
5GFLOPS/core for i7
200s
for 1e5 counts
634y


## results
- the implmentation of above just produces this divergence of the CPD which does not produce a thermal density dist
- seeding with some uncorrelated gaussians does not help

## future thoughts
- i think the answer may lie in the usage og the g2 function
- this algo seems to say that the G^(2)_corr.  is equal to gauss(\delta x,\sigma)+1
- in reality it is the ratio of G^(2)_corr./G^(2)_uncorr. =gauss(\delta x,\sigma)+1
- the G^(2)_uncorr should be related to the  http://mathworld.wolfram.com/NormalDifferenceDistribution.html
- so it is wrong to call the g^(2) the probability of detection relative to the (uncorr) themal distribution


