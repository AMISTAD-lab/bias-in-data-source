import math
from mpmath import *
from data_binarizer import *


def s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs):
    binom_dict = {}
    flipped = False
    bcounts = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(count_vector)
    i = bcounts[0]
    p_b = binary_hypothesis[0]

    def stirl_approx(k):
            return sqrt(n)*mpf(n/math.e)**mpf(n)/(sqrt(2*math.pi*k*(n-k))*mpf(k/math.e)**k*mpf((n-k)/math.e)**(n-k))
    def binom(p):
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
        lo = (s - 10**cur_pow10)
        lo_binom = binom(lo*p_b)
        if lo_binom < alpha:
            return s
        else:
            return buddy(lo)
            
    for x in range(sigfigs):
        cur_val = buddy(cur_val)
        cur_pow10 -= 1

    return (round(cur_val,abs(int(cur_pow10-1))),flipped)


