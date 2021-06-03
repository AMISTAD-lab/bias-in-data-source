import math
import numpy as np
from scipy.optimize import *

def binary_q_finder(hypothesis, p_lowerbound):
    h = np.array([hypothesis[0], hypothesis[1]])
    # Bounds for q[0] and q[1]
    bounds = Bounds([0, 1], [0, 1]) 
    # Since we must have q[0] + q[1] = 1
    linear_constraint = LinearConstraint([1, 1], [0.99], [1.01], keep_feasible=True)

def f(q, h):
    """
    Loss function (using sum of residuals squared)
    which we want to minimize
    """
    return (q[0]-h[0])**2 + (q[1]-h[1])**2

def f_der(q, h):
    """
    Jacobian of f
    """
    return [2*(q[0]-h[0]), 2*(q[1]-h[1])]

def f_hess(q):
    """
    Hessian of f
    """
    return np.array([[2,0], [0,2]])

