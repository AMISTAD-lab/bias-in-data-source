import math
from mg_calculator_event_based import *

def univariate_sc_test(observation, value_list, alpha, hypothesis=[]):
    """
    Conducts a specified complexity hypothesis test 
    for a proposed discrete univariate distribution, using counts.
    Note that our S-method assumes that the uniform distribution is rejected.
    Example: Rolling a weighted die
        observation = 10*[1] + [2,3,4,5,6]
        value_list = [1,2,3,4,5,6]
        hypothesis = [0.1, 0.2, 0.2, 0.4, 0.05, 0.05]
        alpha = 0.05
    Example: Flipping a fair coin
        observation = 20*['H']
        value_list = ['H', 'T']
        hypothesis = [0.5, 0.5]
        alpha = 0.01
    """
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    obs_counts = [observation.count(x) for x in value_list]
    num_bins = len(value_list) #should == len(hyp) == len(obs_counts)
    #now |x| == num of counts 
    norm_scriptx = math.comb(len(observation)+num_bins-1, num_bins-1)
    u = 1/norm_scriptx #probably
    #don't forget the whole-number limitations on this calculator
    mg = mg_calculator_event_based(obs_counts, hyp)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    s_lowerbound = alpha*nu/(r*u)
    p_lowerbound = s_lowerbound*u
    #this should be the PMF, if i've calculated it correctly
    p = (math.factorial(len(observation))/prod([math.factorial(x) for x in obs_counts]))\
    * prod([hyp[x]**obs_counts[x] for x in range(num_bins)]) 
    if p < p_lowerbound:
        reject = True
        print("Proposed distribution rejected at alpha = " + str(alpha) + ". p(x) = " + str(p) + ". s*u(x) = " + str(p_lowerbound) + ".")
    else:
        reject = False
        print("Proposed distribution failed to reject at alpha = " + str(alpha) + ". p(x) = " + str(p) + ". s*u(x) = " + str(p_lowerbound) + ".")
    return (p, s_lowerbound*u, reject)