import math
import numpy as np
from scipy.optimize import *

# 6.575914760628984e-09
def test_binary_q_finder(observation, hypothesis, p_lowerbound):
    h = np.array(hypothesis)
    q0_count = observation.count(0)
    q1_count = observation.count(1)
    def objective(q):
        return (q[0]-h[0])**2 + (q[1]-h[1])**2

    def lin_cons(q):
        return q[0] + q[1] - 1

    def nonlin_cons(q):
        return q[0]**q0_count * q[1]**q1_count - p_lowerbound

    cons1 = ({'type': 'eq', 'fun': lin_cons})
    cons2 = ({'type': 'ineq', 'fun': nonlin_cons})

    for i in np.linspace(0.1, 1, 10):
        x0 = np.array([i, 1-i])
        # print(x0)
        if x0[0]**q0_count * x0[1]**q1_count >= p_lowerbound:
            res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], options={'disp': True})
            return res.x
        else:
            pass
    print("Minimization failed")

def test_q_finder(observation, value_list, hypothesis, p_lowerbound):
    h = np.array(hypothesis)
    value_count = {}
    for value in value_list:
        value_count[value] = observation.count(value)
    def objective(q):
        residuals_squared_list = [(q[i]-h[i])**2 for i in range(len(value_list))]
        print(sum(residuals_squared_list))
        return sum(residuals_squared_list)
    def lin_cons(q):
        return sum(q)-1
    def nonlin_cons(q):
        q_to_count_list = [q[i]**value_count[value_list[i]] for i in range(len(value_list))]
        return math.prod(q_to_count_list) - p_lowerbound
    cons1 = ({'type': 'eq', 'fun': lin_cons})
    cons2 = ({'type': 'ineq', 'fun': nonlin_cons})
    bounds = Bounds(len(value_list)*[0], len(value_list)*[1])
    x0 = np.array(len(value_list)*[1])
    res = minimize(objective, x0, method='SLSQP', constraints=[cons1, cons2], bounds=bounds, options={'ftol': 1e-20,'disp': True})
    return res.x
    