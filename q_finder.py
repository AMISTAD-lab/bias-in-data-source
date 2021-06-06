import math
import numpy as np
import decimal
from scipy.optimize import *
from scipy.special import rel_entr
import warnings

# Used to suppress useless Scipy warnings
# Disable this if you're rigorously testing
warnings.filterwarnings("ignore")

def q_finder_trust_constr(observation, value_list, hypothesis, p_lowerbound):
    # Uses original hypothesis as proposed distribution P
    p = np.array(hypothesis)
    value_count = {}
    for value in value_list:
        value_count[value] = observation.count(value)
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    linear_constraint = LinearConstraint(len(value_list)*[1], [1], [1], keep_feasible=False)
    # Q(x) >= p_lowerbound = s_lowerbound*u(x)
    def nonlin_cons(q):
        q_to_count_list = [q[i]**value_count[value_list[i]] for i in range(len(value_list))]
        return math.prod(q_to_count_list)
    nonlinear_constraint = NonlinearConstraint(nonlin_cons, lb=p_lowerbound, ub=np.inf, keep_feasible=True)
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(value_list)*[0], len(value_list)*[1], keep_feasible=True)
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint, nonlinear_constraint], 
                   options={'gtol': 1e-100, 'xtol': 1e-100, 'maxiter': 1000000, 'verbose': 0}, bounds=bounds)
    return res.x

def q_finder_slsqp(observation, value_list, hypothesis, p_lowerbound):
    # Proposed distribution, p
    p = np.array(hypothesis)
    value_count = {}
    for value in value_list:
        value_count[value] = observation.count(value)
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    def lin_cons(q):
        return sum(q)-1
    # Q(x) >= p_lowerbound = s_lowerbound*u(x)
    def nonlin_cons(q):
        q_to_count_list = [q[i]**value_count[value_list[i]] for i in range(len(value_list))]
        return math.prod(q_to_count_list) - p_lowerbound
    cons1 = ({'type': 'eq', 'fun': lin_cons})
    cons2 = ({'type': 'ineq', 'fun': nonlin_cons})
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(value_list)*[0], len(value_list)*[1])
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], bounds=bounds, 
                   options={'maxiter': 1000000, 'ftol': 1e-100000, 'eps': 1.5e-7, 'disp': False})
    return res.x

# This function overshoots because I don't think
# trust-constr allows for different step sizes
def q_finder_grad_ascent(observation, value_list, hypothesis, p_lowerbound):
    value_count = {}
    for value in value_list:
        value_count[value] = observation.count(value)
    def objective(q):
        q_to_count_list = [q[i]**value_count[value_list[i]] for i in range(len(value_list))]
        q_x = math.prod(q_to_count_list)
        #print(q)
        #print(q_x)
        if q_x >= p_lowerbound:
            return q_x
        return -math.log(q_x, 2)
    linear_constraint = LinearConstraint(len(value_list)*[1], [1], [1], keep_feasible=False)
    bounds = Bounds(len(value_list)*[0], len(value_list)*[1], keep_feasible=True)
    x0 = np.array(hypothesis)
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint], 
                   options={'maxiter': 10000, 'verbose': 0}, bounds=bounds)
    return res.x

# This function also overshoots unfortunately
def q_finder_grad_ascent2(observation, value_list, hypothesis, p_lowerbound):
    value_count = {}
    for value in value_list:
        value_count[value] = observation.count(value)
    def objective(q):
        q_to_count_list = [q[i]**value_count[value_list[i]] for i in range(len(value_list))]
        q_x = math.prod(q_to_count_list)
        #print(q)
        #print(q_x)
        if q_x >= p_lowerbound:
            return q_x
        return -math.log(q_x, 2)
    def lin_cons(q):
        return sum(q)-1
    cons1 = ({'type': 'eq', 'fun': lin_cons})
    bounds = Bounds(len(value_list)*[0.0001], len(value_list)*[1], keep_feasible=True)
    x0 = np.array(hypothesis)
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1], bounds=bounds, 
                   options={'eps': 1.5e-50, 'finite_diff_rel_step': 1e-10, 'disp': True})
    return res.x
    
def q_finder_projection(data,value_list,hypothesis,sux,step_size = 100000000):
    """uses gradient ascent and projection to find q"""
    counts = [decimal.Decimal(data.count(val)) for val in value_list]
    hyp = [decimal.Decimal(hypval) for hypval in hypothesis]
    step_size=decimal.Decimal(step_size)
    while nb_probability_func(hyp,counts) < sux:
        gradient = nb_probability_grad(hyp,counts)
        temphyp = [hyp[i]+gradient[i] for i in range(len(hyp))]
        norm_vector_scalar = (1-sum(temphyp))/len(temphyp)
        hyp = [temphyp[i] + norm_vector_scalar for i in range(len(temphyp))]
    return hyp

def nb_probability_func(hyp,counts):
    """helper for projection"""
    func = 1
    for i in range(len(hyp)):
        func = func*(hyp[i]**counts[i])
    return func

def nb_probability_grad(hyp,counts):
    """helper for projection"""
    grad = []
    for i in range(len(hyp)):
        hypslice = hyp[:i]+hyp[i+1:]
        countsslice = counts[:i]+counts[i+1:]
        grad += [counts[i]*(hyp[i]**(counts[i]-1))*nb_probability_func(hypslice,countsslice)]
    return grad


def q_finder_make_binary(observation, value_list, p_lowerbound):
    """
    If a distribution is rejected as a plausible
    explanation by univariate_sc_test, then this
    function returns a distribution that is the
    closest distribution to the tested hypothesis
    that is also not rejected as a plausible
    explanation.
    """
    assumed_bias = mode(observation)
    assumed_bias_count = observation.count(assumed_bias)
    not_assumed_bias_count = len(observation) - assumed_bias_count
    x = Symbol('x')
    q_list = list(solveset((x**assumed_bias_count) * ((1-x)**not_assumed_bias_count) - p_lowerbound, x, S.Reals))
    q = min([i for i in q_list if i > 0])
    not_q = 1-q
    closest_plausible_dist = []
    for value in value_list:
        if value != assumed_bias:
            closest_plausible_dist.append(not_q/(len(value_list)-1))
        else:
            closest_plausible_dist.append(q)
    return closest_plausible_dist
