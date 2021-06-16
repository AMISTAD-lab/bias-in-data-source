import math
from mpmath import *
import numpy as np
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, Bounds
from scipy.special import rel_entr
import warnings

# Used to suppress useless Scipy warnings
# Disable this if you're rigorously testing
warnings.filterwarnings("ignore")

def q_finder_slsqp(observation, value_list, hypothesis, p_lowerbound):
    # Proposed distribution, p
    p = np.array(hypothesis)
    count_vector = [observation.count(value) for value in value_list]
    # These are for calculating the multinomial PMF
    n_fac = math.factorial(len(observation))
    #print(n_fac)
    prod_x_fac = math.prod([math.factorial(x) for x in count_vector])
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        print(q)
        #print(sum(rel_entr(q, p)))
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    def lin_cons(q):
        return sum(q)-1
    # Q(x) >= p_lowerbound = s_lowerbound*u(x)
    def nonlin_cons(q):
        # This is q0^x0 * q1^x1 * ... * qn^xn
        prod_q = math.prod([q[i]**count_vector[i] for i in range(len(count_vector))])
        #big = n_fac*prod_q/prod_x_fac
        big = log(n_fac) + log(prod_q) - log(prod_x_fac)
        #print(log(prod_q))
        #print(big - log(p_lowerbound))
        return big - log(p_lowerbound)
    cons1 = ({'type': 'eq', 'fun': lin_cons})
    cons2 = ({'type': 'ineq', 'fun': nonlin_cons})
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(count_vector)*[0], len(count_vector)*[1])
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], bounds=bounds, 
                   options={'maxiter': 1000000, 'ftol': 1e-100000, 'eps': 1.5e-7, 'disp': True})
    return res.x

def q_finder_trust_constr(observation, value_list, hypothesis, p_lowerbound):
    # Uses original hypothesis as proposed distribution P
    p = np.array(hypothesis)
    count_vector = [observation.count(value) for value in value_list]
    n_fac = mpf(math.factorial(len(observation)))
    prod_x_fac = math.prod([math.factorial(x) for x in count_vector])
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    linear_constraint = LinearConstraint(len(value_list)*[1], [1], [1], keep_feasible=False)
    # Q(x) >= p_lowerbound = s_lowerbound*u(x)
    def nonlin_cons(q):
        prod_q = math.prod([q[i]**count_vector[i] for i in range(len(count_vector))])
        return n_fac*prod_q/prod_x_fac
    nonlinear_constraint = NonlinearConstraint(nonlin_cons, lb=p_lowerbound, ub=np.inf, keep_feasible=True)
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(value_list)*[0], len(value_list)*[1], keep_feasible=True)
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint, nonlinear_constraint], 
                   options={'gtol': 1e-100, 'xtol': 1e-100, 'maxiter': 1000000, 'verbose': 0}, bounds=bounds)
    return res.x