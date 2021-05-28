from scipy.special import rel_entr
import math


"""I created a version of mgbinary to work in nonbinary situations.
I assume we are given a dict of probability of each element occuring in the 
event, which is prob_dict. If the element is not given a probability, then
its probability is set to its probability of occuring in the given 
distribution. I have not yet figured out how to compute mg based on this, but I 
wanted to make sure the code I have so far makes sense."""
def mgnonbinary(data,prob_dict):
    numofvals = len(data)
    if sum(prob_dict.values()) != 1:
        raise Exception("Given Probability Distribution does not sum to 1")
	for x in data:
		if x not in prob_dict.keys():
			probofx = data.count(x)/numofvals
            prob_dict[x] = probofx
        
    uni_dist = len(prob_dict.keys())*[1/len(prob_dict.keys())]
	obs_dist = [data.count(i)/numofvals for i in prob_dict.keys()]
    
    min_kl = sum(rel_entr(uni_dist, obs_dist))
    mg = 0
	
	