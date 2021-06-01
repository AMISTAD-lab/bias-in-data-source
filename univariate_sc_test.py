import math
from itertools import *
from scipy.special import rel_entr
from sympy.solvers import solve
from sympy import Symbol
from statistics import *

def closest_plausible_explanation(p_lowerbound, assumed_bias_count, not_assumed_bias_count):
    q = Symbol('q', real=True, positive=True)
    solutions = solve((q**assumed_bias_count) * ((1-q)**not_assumed_bias_count) - p_lowerbound, q)
    return solutions

def univariate_sc_test(observation, value_list, hypothesis, alpha):
    """
    Conducts a specified complexity hypothesis test 
    for a proposed discrete univariate distribution
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
    norm_scriptx = len(value_list)**len(observation)
    u = 1/norm_scriptx
    nu = norm_scriptx/mg(observation, value_list)
    r = norm_scriptx*(1+math.log(norm_scriptx))
    s_lowerbound = alpha*nu/(r*u)
    p_lowerbound = s_lowerbound*u
    p = math.prod([hypothesis[value_list.index(i)] for i in observation])
    if p < p_lowerbound:
        reject = True
        print("Proposed distribution rejected at alpha = " + str(alpha) + ". p(x) = " + str(p) + ". s*u(x) = " + str(p_lowerbound) + ".")
        assumed_bias = mode(observation)
        assumed_bias_count = observation.count(assumed_bias)
        not_assumed_bias_count = len(observation) - assumed_bias_count
        q = closest_plausible_explanation(p_lowerbound, assumed_bias_count, not_assumed_bias_count)[0]
        not_q = (1-q)/(len(value_list)-1)
        closest_plausible_dist = []
        for value in value_list:
            if value != assumed_bias:
                closest_plausible_dist.append(not_q)
            else:
                closest_plausible_dist.append(q)
        print("Closest plausible distribution: " + str(closest_plausible_dist))
    else:
        reject = False
        print("Proposed distribution failed to reject at alpha = " + str(alpha) + ". p(x) = " + str(p) + ". s*u(x) = " + str(p_lowerbound) + ".")
    return (p, s_lowerbound*u, reject)

def uniform_dist_sc_test(observation, value_list, alpha):
    """
    Conducts a specified complexity hypothesis test 
    for a uniform discrete univariate distribution
    Example: Flippin a fair coin
        observation: 20*['H']
        value_list = ['H', 'T']
        alpha = 0.05
    """
    norm_scriptx = len(value_list)**len(observation)
    u = 1/norm_scriptx
    nu = norm_scriptx/mg(observation, value_list)
    r = norm_scriptx*(1+math.log(norm_scriptx))
    kardis = r*u/nu
    s_lowerbound = alpha*nu/(r*u)
    if kardis < alpha:
        reject = True
        print("Uniform distribution rejected at alpha = "+ str(alpha) +". Kardis = " + str(kardis)+ ". s >", s_lowerbound)
    else:
        reject = False
        print("Uniform distribution not rejected at alpha = "+ str(alpha) +". Kardis = " + str(kardis)+ ". s >", s_lowerbound)
    return (kardis, s_lowerbound, reject)

def mg(observation, value_list):
    """
    Calculates Mg(x) for a given observation and list of 
    possible values for a given discrete random variable
    utilizing KL Divergence as a difference measure between
    two distributions
    """
    freq_dict = {}
    for value in value_list:
        freq_dict[value] = observation.count(value)
    uni_dist = len(freq_dict.keys())*[1/len(freq_dict.keys())]
    obs_dist = [observation.count(i)/len(observation) for i in freq_dict.keys()]
    min_kl = sum(rel_entr(uni_dist, obs_dist))
    mg = 0
    scriptx = scriptx_generator_helper(scriptx_generator(len(value_list), len(observation)), len(value_list))
    for event in scriptx:
        test_dist = [event[i]/len(observation) for i in range(len(value_list))]
        if sum(rel_entr(uni_dist, test_dist)) >= min_kl:
            mg += math.factorial(len(observation)) // math.prod(list(map(lambda x: math.factorial(x), event)))
    return mg

def scriptx_generator(num_vals, observation_length, current_vals=[]):
    """
    Generates a representation of scriptx
    for a discrete random variable
    Notes:
        This representation of scriptx is based
        "bins" of values, where scriptx is a list
        of lists. Each list in scriptx is a
        possible frequency distribution. However,
        this function itself does not do the grouping
        of these distributions into lists. Rather, this
        function outputs just a single list of numbers
        that are grouped into distributions by the helper
        function below.
    Example: Rolling a fair die 6 times
        Some of the possible frequency distributions
        in scriptx are [1,1,1,1,1,1], [0,6,0,0,0,0],
        and [1,0,3,2,0,0]. These would show up in scriptx
        as the single list 
        [1,1,1,1,1,1,0,6,0,0,0,0,1,0,3,2,0,0]
        and would be grouped into nice lists by the helper
        function below.
    """
    if observation_length == 0:
        return current_vals + [0]*num_vals
    elif num_vals == 1:
        return current_vals + [observation_length]
    else:
        scriptx = []
        for i in range(observation_length+1):
            scriptx += current_vals + scriptx_generator(num_vals-1, observation_length-i, current_vals=[i])
        return scriptx

def scriptx_generator_helper(ungrouped_perm_list, num_vals):
    """
    Helper function for scriptx_generator
    Notes:
        scriptx generator only generates a single
        list of numbers which are in the correct order,
        but still need to be grouped into frequency
        distributions in order to be interpretable.
        This function does this grouping.
    """
    grouped_perm_list = []
    for i in range(0, len(ungrouped_perm_list), num_vals):
        grouped_perm_list.append(ungrouped_perm_list[i:i+num_vals])
    return grouped_perm_list
