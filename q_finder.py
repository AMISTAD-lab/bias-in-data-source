import math
import numpy as np
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, Bounds
from scipy.special import rel_entr
import warnings

# Used to suppress useless Scipy warnings
# Disable this if you're rigorously testing
warnings.filterwarnings("ignore")

def q_finder_main_slsqp(count_vector, hypothesis, p_lowerbound):
    """
    Returns the closest plausible distribution to the 'hypothesis' distribution
    using sequential quadratic programing. 
    'count_vector' contains the quantity of each value in the data. 
    'p_lowerbound' is the lower bound on plausible probabilities - any
    distribution which produces a lower sequence probability
    may be instantly rejected.
    """
    # Proposed distribution, p
    p = np.array(hypothesis)
    # These are for calculating the multinomial PMF
    log_n_fac = math.log(math.factorial(sum(count_vector)))
    log_prod_x_fac = math.log(math.prod([math.factorial(x) for x in count_vector]))
    constraint_constant = math.log(p_lowerbound) - log_n_fac + log_prod_x_fac
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    def lin_cons(q):
        return sum(q)-1
    # Q(x) >= p_lowerbound = s_lowerbound*u(x)
    def nonlin_cons(q):
        log_prod_q = sum([count_vector[i]*math.log(q[i]) for i in range(len(count_vector))])
        constraint = log_prod_q - constraint_constant
        return constraint
    cons1 = ({'type': 'eq', 'fun': lin_cons})
    cons2 = ({'type': 'ineq', 'fun': nonlin_cons})
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(count_vector)*[1e-10], len(count_vector)*[1])
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(count_vector)*[1])
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], bounds=bounds, 
                   options={'disp': False})
    return res.x

def q_finder_main_trust_constr(count_vector, hypothesis, p_lowerbound):
    """
    Returns the closest plausible distribution to the 'hypothesis' distribution
    using trust constraint gradient descent. 
    'count_vector' contains the quantity of each value in the data. 
    'p_lowerbound' is the lower bound on plausible probabilities - any
    distribution which produces a lower sequence probability
    may be instantly rejected.
    """
    p = np.array(hypothesis)
    log_n_fac = math.log(math.factorial(sum(count_vector)))
    log_prod_x_fac = math.log(math.prod([math.factorial(x) for x in count_vector]))
    constraint_constant = math.log(p_lowerbound) - log_n_fac + log_prod_x_fac
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    linear_constraint = LinearConstraint(len(count_vector)*[1], [1], [1], keep_feasible=False)
    # log(Q(x)) - constraint_constant >= 0
    def nonlin_cons(q):
        log_prod_q = sum([count_vector[i]*math.log(q[i]) for i in range(len(count_vector))])
        constraint = log_prod_q - constraint_constant
        return constraint
    nonlinear_constraint = NonlinearConstraint(nonlin_cons, lb=0, ub=np.inf, keep_feasible=True)
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(count_vector)*[0], len(count_vector)*[1], keep_feasible=True)
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(count_vector)*[0.9])
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint, nonlinear_constraint], 
                   options={'gtol': 1e-100, 'xtol': 1e-100, 'maxiter': 1000000, 'verbose': 0}, bounds=bounds)
    return res.x

def binary_q_finder_main(binary_count_vector, hypothesis, alpha):
    """
    Acts similarly to 'q_finder_main', but is strictly for binary situations.
    'hypothesis' is the originally proposed distribution, 
    and 'alpha' is the given alpha level. 
    """
    p = np.array(hypothesis)
    n = sum(binary_count_vector)
    k = binary_count_vector[0]
    constraint_constant = math.log(alpha) - math.log(n+1) - math.log(math.comb(n,k))
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
