import math
from mg_calculator_count_based import *
from mpmath import *
from q_finder_count_based import *
from observation_count import *

def uniform_dist_kardis_test(observation, value_list, alpha):
    obs_counts = observation_count(observation,value_list)
    num_bins = len(value_list) #should == len(hyp) == len(obs_counts)
    #now |x| == num of counts 
    norm_scriptx = math.comb(len(observation)+num_bins-1, num_bins-1)
    u_hyp = [1/num_bins]*num_bins
    u = mpf(math.factorial(len(observation))) / math.prod([mpf(math.factorial(x)) for x in obs_counts])\
    * math.prod([mpf(u_hyp[x])**mpf(obs_counts[x]) for x in range(num_bins)])
    mg = mg_calculator_event_based(obs_counts, u_hyp)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    kardis = r*u/nu
    if kardis < alpha:
        reject = True
        print("Uniform distribution rejected at alpha = "+ str(alpha) + ". p(x) = " + str(u) + ". Kardis = " + str(kardis))
    else:
        reject = False
        print("Uniform distribution not rejected at alpha = "+ str(alpha)  + ". p(x) = " + str(u) +". Kardis = " + str(kardis))
    return (kardis, reject)

def univariate_kardis_test(observation, value_list, alpha, hypothesis=[]):
    #resort to this test if uniform dist is not rejected?
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    obs_counts = observation_count(observation,value_list)
    num_bins = len(value_list) #should == len(hyp) == len(obs_counts)
    #now |x| == num of counts 
    norm_scriptx = math.comb(len(observation)+num_bins-1, num_bins-1)
    #don't forget the whole-number limitations on this calculator
    mg = mg_calculator_event_based(obs_counts, hyp)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    p = mpf(math.factorial(len(observation))) / math.prod([mpf(math.factorial(x)) for x in obs_counts])\
    * math.prod([mpf(hyp[x])**mpf(obs_counts[x]) for x in range(num_bins)])
    #print(p)
    kardis = r*p/nu
    if kardis < alpha:
        reject = True
        print("Proposed distribution rejected at alpha = "+ str(alpha)  + ". p(x) = " + str(p) +". Kardis = " + str(kardis))
    else:
        reject = False
        print("Proposed distribution not rejected at alpha = "+ str(alpha)  + ". p(x) = " + str(p) +". Kardis = " + str(kardis))
    return (kardis, reject)
