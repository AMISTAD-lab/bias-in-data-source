import math
import numpy as np
from scipy.optimize import *
from scipy.special import rel_entr

def binary_q_finder(observation, hypothesis, p_lowerbound):
    h = np.array([hypothesis[0], hypothesis[1]])
    # q0 is heads, q1 is tails, sorry for the confusion
    # for some reason, the program doesn't work if you
    # flip them :(
    q0_count = observation.count(1)
    q1_count = observation.count(0)
    def f(q):
        """
        Loss function (using sum of residuals squared)
        which we want to minimize
        """
        return (q[0]-h[0])**2 + (q[1]-h[1])**2
    def f_der(q):
        """
        Jacobian of f
        """
        return [2*(q[0]-h[0]), 2*(q[1]-h[1])]
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
        q[0]^x[0] * q[1]^x[1] >= p_lowerbound
        """
        return (q[0]**q0_count) * (q[1]**q1_count)
    def cons_J(q):
        """
        Jacobian of the constraint function
        """
        if q0_count == 0:
            return [0, q1_count*(q[1]**(q1_count-1))]
        elif q1_count == 0:
            return [q0_count*(q[0]**(q0_count-1)), 0]
        else:
            return [q0_count*(q[0]**(q0_count-1))*(q[1]**q1_count), q1_count*(q[0]**q0_count)*(q[1]**(q1_count-1))]
    def cons_H(q, v):
        """
        Hessian of the constraint function
        """
        if q0_count == 0:
            return v[0]*np.array([[0,0], [0, q1_count*(q1_count-1)*(q[1]**(q1_count-2))]])
        elif q1_count == 0:
            return v[0]*np.array([[q0_count*(q0_count-1)*(q[0]**(q0_count-2)),0], [0,0]])
        else:
            return v[0]*np.array([[(q0_count*(q0_count-1))*(q[0]**(q0_count-2))*(q[1]**q1_count), (q0_count*q1_count)*(q[0]**(q0_count-1))*(q[1]**(q1_count-1))],
                                  [(q0_count*q1_count)*(q[0]**(q0_count-1))*(q[1]**(q1_count-1)), (q1_count*(q1_count-1))*(q[0]**q0_count)*(q[1]**(q1_count-2))]])
    # Constructs NonLinearConstraint object
    nonlinear_constraint = NonlinearConstraint(cons_f, lb=p_lowerbound, ub=np.inf, jac=cons_J, hess=cons_H, keep_feasible=True)
    
    # This is supposed to minimize f subject to the constraints
    # We need to loop through different x0 values because the minimize
    # function requires an x0 that meets the constraints
    for i in np.linspace(0.1, 1, 10):
        x0 = np.array([i, 1-i])
        # print(x0)
        try:
            res = minimize(f, x0, method='trust-constr', jac=f_der, hess=f_hess, 
                   constraints=[linear_constraint, nonlinear_constraint], options={'verbose': 2}, bounds=bounds)
            return res.x
        except:
            pass
    return "Minimization failed"

def binary_q_finder_kl(observation, hypothesis, p_lowerbound):
    h = np.array([hypothesis[0], hypothesis[1]])
    # q0 is heads, q1 is tails, sorry for the confusion
    # for some reason, the program doesn't work if you
    # flip them :(
    q0_count = observation.count(0)
    q1_count = observation.count(1)
    def f(q):
        """
        Loss function (using sum of residuals squared)
        which we want to minimize
        """
        return sum(rel_entr(q, h))
    def f_der(q):
        """
        Jacobian of f
        """
        return [2*(q[0]-h[0]), 2*(q[1]-h[1])]
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
        q[0]^x[0] * q[1]^x[1] >= p_lowerbound
        """
        return (q[0]**q0_count) * (q[1]**q1_count)
    def cons_J(q):
        """
        Jacobian of the constraint function
        """
        if q0_count == 0:
            return [0, q1_count*(q[1]**(q1_count-1))]
        elif q1_count == 0:
            return [q0_count*(q[0]**(q0_count-1)), 0]
        else:
            return [q0_count*(q[0]**(q0_count-1))*(q[1]**q1_count), q1_count*(q[0]**q0_count)*(q[1]**(q1_count-1))]
    def cons_H(q, v):
        """
        Hessian of the constraint function
        """
        if q0_count == 0:
            return v[0]*np.array([[0,0], [0, q1_count*(q1_count-1)*(q[1]**(q1_count-2))]])
        elif q1_count == 0:
            return v[0]*np.array([[q0_count*(q0_count-1)*(q[0]**(q0_count-2)),0], [0,0]])
        else:
            return v[0]*np.array([[(q0_count*(q0_count-1))*(q[0]**(q0_count-2))*(q[1]**q1_count), (q0_count*q1_count)*(q[0]**(q0_count-1))*(q[1]**(q1_count-1))],
                                  [(q0_count*q1_count)*(q[0]**(q0_count-1))*(q[1]**(q1_count-1)), (q1_count*(q1_count-1))*(q[0]**q0_count)*(q[1]**(q1_count-2))]])
    # Constructs NonLinearConstraint object
    nonlinear_constraint = NonlinearConstraint(cons_f, lb=p_lowerbound, ub=np.inf, jac=cons_J, hess=cons_H, keep_feasible=True)
    
    # This is supposed to minimize f subject to the constraints
    # We need to loop through different x0 values because the minimize
    # function requires an x0 that meets the constraints
    for i in np.linspace(0.1, 1, 10):
        x0 = np.array([i, 1-i])
        print(x0)
        try:
            res = minimize(f, x0, method='trust-constr', 
                   constraints=[linear_constraint, nonlinear_constraint], options={'verbose': 2}, bounds=bounds)
            return res.x
        except:
            pass
    return "Minimization failed"