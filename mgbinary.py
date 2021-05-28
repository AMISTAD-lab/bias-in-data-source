import math
from scipy.special import rel_entr

# This comment is to help understand how the rel_entr function works
"""
# define distributions
p = [0.5, 0.5]
q = [1, 0]
r = [0.8, 0.2]
# calculate (P || Q)
kl_pq = rel_entr(p, q)
print('KL(P || Q): %.3f nats' % sum(kl_pq))
# calculate (Q || P)
kl_pr = rel_entr(p, r)
print('KL(P || R): %.3f nats' % sum(kl_pr))
"""

def mgbinary(data):
    uni_dist = [0.5, 0.5]
    obs_dist = [data.count(1)/len(data), data.count(0)/len(data)]
    min_kl = sum(rel_entr(uni_dist, obs_dist))
    mg = 0
    for i in range(len(data)+1):
        test_obs = (i)*[1] + (len(data)-i)*[0]
        test_dist = [test_obs.count(1)/len(test_obs), test_obs.count(0)/len(test_obs)]
        if sum(rel_entr(uni_dist, test_dist)) >= min_kl:
            #print(test_obs)
            mg += math.comb(len(test_obs), test_obs.count(1))
    return mg

def mg(data, numvalues):
    freq_dict = {}
    for entry in data:
        if entry in freq_dict:
            freq_dict[entry] += 1
        else:
            freq_dict[entry] = 1
    uni_dist = len(freq_dict.keys())*[1/len(freq_dict.keys())]
    obs_dist = [data.count(i)/len(data) for i in freq_dict.keys()]
    min_kl = sum(rel_entr(uni_dist, obs_dist))
    mg = 0
    


