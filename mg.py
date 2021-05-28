import math
from scipy.special import rel_entr
from itertools import *

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

def mg(data, value_list):
    freq_dict = {}
    """
    for entry in data:
        if entry in freq_dict:
            freq_dict[entry] += 1
        else:
            freq_dict[entry] = 1
    """
    for value in value_list:
        freq_dict[value] = data.count(value)

    uni_dist = len(freq_dict.keys())*[1/len(freq_dict.keys())]
    obs_dist = [data.count(i)/len(data) for i in freq_dict.keys()]
    min_kl = sum(rel_entr(uni_dist, obs_dist))
    mg = 0
    for perm in all_perms(value_list, len(data)):
        #print(perm)
        test_dist = [perm.count(i)/len(data) for i in value_list]
        if sum(rel_entr(uni_dist, test_dist)) >= min_kl:
            mg += 1
    return mg

# Helper function for mg
def all_perms(value_list, k):
    full_blocks = k//len(value_list)
    val_dict = {}
    for i in range(full_blocks):
        perm_list = list(product(value_list, repeat=len(value_list)))
        for j in range(len(perm_list)):
            perm_list[j] = list(perm_list[j])
        val_dict[i] = perm_list
    if k%len(value_list) != 0:
        perm_list = list(product(value_list, repeat=k%len(value_list)))
        for i in range(len(perm_list)):
            perm_list[i] = list(perm_list[i])
        val_dict[full_blocks] = perm_list
    block_list = []
    for key in val_dict.keys():
        block_list.append(val_dict[key])
    total_perm_list = list(product(*block_list))
    for i in range(len(total_perm_list)):
        total_perm_list[i] = list(chain.from_iterable(total_perm_list[i]))
    return total_perm_list





