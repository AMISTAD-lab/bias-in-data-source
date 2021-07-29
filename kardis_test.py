import math
from mg_calculator import *
from mpmath import *

def kardis_test(counts, alpha, hypothesis):
    """
    Performs a multinomial kardis test upon the provided count vector 'counts'.
    'alpha' is the chosen alpha level, and 'hypothesis' is a list of 
    proposed probabilities for each value, 
    where counts[i] corresponds to hypothesis[i].
    Returns the kardis, a boolean corresponding to rejection, 
    as well as three values purely for use in the main function that calls this.
    """
    num_bins = len(counts)
    n = sum(counts)
    norm_scriptx = math.comb(n+num_bins-1, num_bins-1)
    mg = mg_calculator(counts, hypothesis)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    h = mpf(math.factorial(n)) / math.prod([mpf(math.factorial(x)) for x in counts])\
    * math.prod([mpf(hypothesis[x])**mpf(counts[x]) for x in range(num_bins)])
    kardis = r*h/nu
    if kardis < alpha:
        reject = True
    else:
        reject = False
    return (kardis, reject, r, nu, h)