import math
import numpy as np
import decimal
from scipy.optimize import *

def q_finder_slsqp(observation, value_list, hypothesis, p_lowerbound):
    # Proposed distribution, p
    p = np.array(hypothesis)
    value_count = {}
    for value in value_list:
        value_count[value] = observation.count(value)
    def objective(q):
        """
        Objective function which we want to minimize
        Defined by sum of square residuals between q and p
        """
        residuals_squared_list = [(q[i]-p[i])**2 for i in range(len(value_list))]
        return sum(residuals_squared_list)
    def lin_cons(q):
        """
        Constraint that the conditional probabilities
        of Q must sum to 1
        """
        return sum(q)-1
    def nonlin_cons(q):
        """
        Constraint that Q(x) >= s*u(x)
        """
        q_to_count_list = [q[i]**value_count[value_list[i]] for i in range(len(value_list))]
        return math.prod(q_to_count_list) - p_lowerbound
    cons1 = ({'type': 'eq', 'fun': lin_cons})
    cons2 = ({'type': 'ineq', 'fun': nonlin_cons})
    # For each conditional probability of Q, q[i]
    # we must have 0 <= q[i] <= 1
    bounds = Bounds(len(value_list)*[0], len(value_list)*[1])
    # Bogus initial guess of all 1's
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], bounds=bounds, 
                   options={'maxiter': 1000000, 'ftol': 1e-1000, 'eps': 1.5e-7, 'disp': False})
    return res.x

def q_finder_trust_constr(observation, value_list, hypothesis, p_lowerbound):
    p = np.array(hypothesis)
    value_count = {}
    for value in value_list:
        value_count[value] = observation.count(value)
    def objective(q):
        residuals_squared_list = [(q[i]-p[i])**2 for i in range(len(value_list))]
        return sum(residuals_squared_list)
    linear_constraint = LinearConstraint(len(value_list)*[1], [1], [1], keep_feasible=False)
    def nonlin_cons(q):
        q_to_count_list = [q[i]**value_count[value_list[i]] for i in range(len(value_list))]
        return math.prod(q_to_count_list)
    nonlinear_constraint = NonlinearConstraint(nonlin_cons, lb=p_lowerbound, ub=np.inf, keep_feasible=True)
    bounds = Bounds(len(value_list)*[0], len(value_list)*[1], keep_feasible=True)
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint, nonlinear_constraint], 
                   options={'gtol': 1e-100, 'xtol': 1e-100, 'maxiter': 1000000, 'verbose': 1}, bounds=bounds)
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
