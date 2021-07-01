import math
from mpmath import *
from data_binarizer import *


def s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs):
    """
    For the tightest-bound binary function. Finds the coefficient necessary 
    on the selected value probability to produce a plausible distribution. 
    'count_vector' is the original non-binary count representation of the data,
    'value_list' contains all the values a data point may take, 
    'selected_value_list' is the collection of values the user believes
    the data is biased towards. 'alpha' is the requested alpha-level, 
    'binary_hypothesis' contains the probability of producing a
    selected value or a not-selected value,
    and 'sigfigs' is the requested number of significant figures.
    """
    binom_dict = {}
    flipped = False
    bcounts = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(count_vector)
    i = bcounts[0]
    p_b = binary_hypothesis[0]

    def stirl_approx(k):
            """
            Approximates nCk using Stirling's approximation.
            """
            return sqrt(n)*mpf(n)**mpf(n)/(sqrt(2*math.pi*k*(n-k))*mpf(k)**k*mpf(n-k)**(n-k))
    def binom(p):
        """
        Returns the binomial one-sided (greater) tail probability.
        """
        if p in binom_dict.keys():
            return binom_dict[p]
        else:
            prob_more_extreme = sum([stirl_approx(k)*(mpf(p)**k)*(mpf(1-p)**(n-k)) for k in range(i,n)])
            prob_more_extreme += mpf(p)**n
            binom_dict[p] = prob_more_extreme
            return prob_more_extreme

    if binom(p_b) >= alpha:
        if binom(p_b) > 1-alpha:
            flipped = True
            i = bcounts[1]
            p_b = binary_hypothesis[1]
        else:
            return 1, flipped
    cur_pow10 = 0
    cur_val = 1/p_b

    def buddy(s):
        """
        Recursively optimizes s' for a single power of 10.
        """
        lo = (s - 10**cur_pow10)
        lo_binom = binom(lo*p_b)
        if lo_binom < alpha:
            return s
        else:
            return buddy(lo)
            
    for x in range(sigfigs):
        cur_val = buddy(cur_val)
        cur_pow10 -= 1

    return (round(cur_val,abs(int(cur_pow10+1))),flipped)


