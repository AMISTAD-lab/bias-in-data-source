import math
import numpy as np
from scipy.optimize import *

"""
HOW TO USE:
    1. Run this file
    2. Type into terminal print(res.x)

See more detailed instructions at the section of this page 
titled "Unconstrained minimization of multivariate scalar functions"
https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html#sequential-least-squares-programming-slsqp-algorithm-method-slsqp

This is for the example of flipping a coin and getting the sequence
20*[1]+5*[0] where 1 is heads and 0 is tails.

Running this through univariate_sc_test w/ hypothesis=[0.5,0.5]
and alpha=0.05 gets s*u(x) = 6.575914760628984e-09.

I at least got this to run, but the optimizer is not obeying
my nonlinear constraint of q[0]^25 * q[1]^5 >= 6.575914760628984e-09 
:(
"""
def f(q):
    """
    Loss function (using sum of residuals squared)
    which we want to minimize
    """
    res = (q[0]-0.5)**2 + (q[1]-0.5)**2
    print('{0:4f}   {1: 3.6f}   {2: 3.6f}   '.format(q[0], q[1], res))
    return res
def f_der(q):
    """
    Jacobian of f
    """
    return [2*q[0]-1, 2*q[1]-1]
def f_hess(q):
    """
    Hessian of f
    """
    return np.array([[2,0], [0,2]])
# Bounds for q[0] and q[1]
bounds = Bounds([0, 1], [0, 1]) 
# Since we must have q[0] + q[1] = 1
linear_constraint = LinearConstraint([1, 1], [0.99], [1.01], keep_feasible=True)
def cons_f(q):
    """
    Constraint function: 
    q[0]^25 * q[1]^5 >= s*u(x) = 6.575914760628984e-09
    """
    return (q[0]**5) * (q[1]**25)
def cons_J(q):
    """
    Jacobian of the constraint function
    """
    return [5*(q[0]**4)*(q[1]**25), 25*(q[0]**5)*(q[1]**24)]
def cons_H(q, v):
    """
    Hessian of the constraint function
    """
    return v[0]*np.array([[20*(q[0]**3)*(q[1]**25), 125*(q[0]**4)*(q[1]**24)],
                            [125*(q[0]**4)*(q[1]**24), 600*(q[0]**5)*(q[1]**23)]])
# Constructs NonLinearConstraint object
nonlinear_constraint = NonlinearConstraint(cons_f, lb=6.575914760628984e-09, ub=np.inf, jac=cons_J, hess=cons_H, keep_feasible=False)
# Starting point for the optimization
x0 = np.array([0.4,0.6])

Nfeval = 1

def callbackF(qi):
    global Nfeval
    print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   '.format(Nfeval, qi[0], qi[1], f(qi)))
    Nfeval += 1

# This is supposed to minimize f subject to the constraints
res = minimize(f, x0, method='trust-constr', jac=f_der, hess=f_hess, 
            constraints=[linear_constraint, nonlinear_constraint], options={'verbose': 0, 'disp': True}, bounds=bounds)
print(res.x)
