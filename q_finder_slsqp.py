import math
import numpy as np
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
        #print(q)
        #print(math.prod([q[i]**value_count[value_list[i]] for i in range(len(value_list))]))
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
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], bounds=bounds, options={'maxiter': 1000000, 'ftol': 1e-100, 'eps': 1.5e-7, 'disp': True})
    # print(math.prod([res.x[i]**value_count[value_list[i]] for i in range(len(value_list))]))
    return res.x