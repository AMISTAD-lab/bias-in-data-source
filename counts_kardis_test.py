import math
from mg_calculator_count_based import *
from mpmath import *

def kardis_test_main(counts, alpha, hypothesis):
    num_bins = len(counts)
    n = sum(counts)
    norm_scriptx = math.comb(n+num_bins-1, num_bins-1)
    mg = mg_calculator(counts, hypothesis)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    h = mpf(math.factorial(n)) / math.prod([mpf(math.factorial(x)) for x in counts])\
    * math.prod([mpf(hypothesis[x])**mpf(counts[x]) for x in range(num_bins)])
    kardis = r*h/nu
    if kardis < alpha:
        reject = True
    else:
        reject = False
    return (kardis, reject, r, nu, h)

def binary_kardis_test_main(counts, alpha, hypothesis):
    n = sum(counts)
    k = counts[0]
    r = n+1
    p = mpf(hypothesis[0])**counts[0] * mpf(hypothesis[1])**counts[1]
    nu = mpf(math.comb(n, k))**(-1)
    kardis = r*p/nu
    if kardis < alpha:
        reject = True
    else:
        reject = False
    return (kardis, reject, r, nu, p)

def uniform_dist_kardis_test(observation, value_list, alpha):
    obs_counts = [observation.count(x) for x in value_list]
    num_bins = len(value_list) #should == len(hyp) == len(obs_counts)
    #now |x| == num of counts 
    norm_scriptx = math.comb(len(observation)+num_bins-1, num_bins-1)
    u_hyp = [1/num_bins]*num_bins
    u = mpf(math.factorial(len(observation))) / math.prod([mpf(math.factorial(x)) for x in obs_counts])\
    * math.prod([mpf(u_hyp[x])**mpf(obs_counts[x]) for x in range(num_bins)])
    mg = mg_calculator(obs_counts, u_hyp)
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
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    obs_counts = [observation.count(x) for x in value_list]
    num_bins = len(value_list) #should == len(hyp) == len(obs_counts)
    #now |x| == num of counts 
    norm_scriptx = math.comb(len(observation)+num_bins-1, num_bins-1)
    #don't forget the whole-number limitations on this calculator
    mg = mg_calculator(obs_counts, hyp)
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