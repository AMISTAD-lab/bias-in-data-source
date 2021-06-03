import math
import numpy as np
from scipy.optimize import *

def binary_q_finder(hypothesis, p_lowerbound):
    h = np.array([hypothesis[0], hypothesis[1]])

def f(q, h):
    return (q[0]-h[0])**2 + (q[1]-h[1])**2

