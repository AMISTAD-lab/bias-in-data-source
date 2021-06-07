import math
from mg import *
import numpy as np
from q_finder import *
from scipy.optimize import *
from scipy.special import rel_entr
from statistics import *
from sympy import Symbol, solveset, S

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
    mg = mg_calculator(observation, value_list)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    s_lowerbound = alpha*nu/(r*u)
    p_lowerbound = s_lowerbound*u
    p = math.prod([hypothesis[value_list.index(i)] for i in observation])
    if p < p_lowerbound:
        reject = True
        print("Proposed distribution rejected at alpha = " + str(alpha) + ". p(x) = " + str(p) + ". s*u(x) = " + str(p_lowerbound) + ".")
        Q = list(q_finder_slsqp(observation, value_list, hypothesis, p_lowerbound))
        print("Closest plausible distribution: " + str(Q))
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
    mg = mg_calculator(observation, value_list)
    nu = norm_scriptx/mg
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