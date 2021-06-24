import math
from mpmath import *
import numpy as np
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, Bounds
from scipy.special import rel_entr
import warnings

# Used to suppress useless Scipy warnings
# Disable this if you're rigorously testing
warnings.filterwarnings("ignore")
# I don't believe this one works just yet
def q_finder_slsqp(observation, value_list, hypothesis, p_lowerbound):
    # Proposed distribution, p
    p = np.array(hypothesis)
    count_vector = [observation.count(x) for x in value_list]
    # These are for calculating the multinomial PMF
    n_fac = mpf(math.factorial(len(observation)))
    prod_x_fac = math.prod([math.factorial(x) for x in count_vector])
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        #print(q)
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    def lin_cons(q):
        return sum(q)-1
    # Q(x) >= p_lowerbound = s_lowerbound*u(x)
    def nonlin_cons(q):
        # This is q0^x0 * q1^x1 * ... * qn^xn
        prod_q = math.prod([mpf(q[i])**mpf(count_vector[i]) for i in range(len(count_vector))])
        try:
            log_prod_q = math.log(prod_q)
        except:
            log_prod_q = -10000
        if math.log(n_fac) == inf:
            log_n_fac = 10000
        else:
            log_n_fac = math.log(n_fac)
        #print(str(log_n_fac)+ ' ' +str(log_prod_q)+ ' ' +str(math.log(prod_x_fac)) + ' ' + str(math.log(p_lowerbound)))
        #return log_n_fac+log_prod_q-math.log(prod_x_fac) - math.log(p_lowerbound)
        return n_fac*prod_q/prod_x_fac - p_lowerbound
    cons1 = ({'type': 'eq', 'fun': lin_cons})
    cons2 = ({'type': 'ineq', 'fun': nonlin_cons})
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(count_vector)*[0], len(count_vector)*[1])
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], bounds=bounds, 
                   options={'disp': True})
    return res.x