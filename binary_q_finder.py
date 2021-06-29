import math
from mpmath import *
import numpy as np
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, Bounds
from scipy.special import rel_entr
import warnings

# Used to suppress useless Scipy warnings
# Disable this if you're rigorously testing
warnings.filterwarnings("ignore")

def binary_q_finder(binary_count_vector, hypothesis, alpha):
    p = np.array(hypothesis)
    n = sum(binary_count_vector)
    k = binary_count_vector[0]
    constraint_constant = math.log(alpha * (n+1)**(-1) * math.comb(n, k)**(-1))
    def objective(q):
        return sum(rel_entr(q, p))
    linear_constraint = LinearConstraint([1,1], [1], [1], keep_feasible=False)
    def nonlin_cons(q):
        return k*math.log(q[0]) + (n-k)*math.log(q[1]) - constraint_constant
    nonlinear_constraint = NonlinearConstraint(nonlin_cons, lb=0, ub=np.inf, keep_feasible=True)
    bounds = Bounds([0,0], [1,1], keep_feasible=True)
    x0 = np.array([0.99,0.99])
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint, nonlinear_constraint], 
                   options={'gtol': 1e-100, 'xtol': 1e-100, 'maxiter': 1000000, 'verbose': 0}, bounds=bounds)
    return res.x