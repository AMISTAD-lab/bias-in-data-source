import math
from mpmath import *
from data_binarizer import *


def s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs):
    bcounts = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(count_vector)
    i = bcounts[0]
    p_b = binary_hypothesis[0]
    def stirl_approx(k):
            return (mpf(n/math.e)**n)/(sqrt(2*math.pi*k*(n-k)) * mpf(k/math.e)**k * mpf((n-k)/math.e)**(n-k))
    def binom(p):
        prob_more_extreme = sum([stirl_approx(k)*(mpf(p)**k)*(mpf(1-p)**(n-k)) for k in range(i,n)])
        prob_more_extreme += mpf(p)**n
        return prob_more_extreme
    if binom(p_b) >= alpha:
        print("fine")
        return 1
    else:
        cur_pow10 = 0
        cur_val = 1
        upper_bound = 1/p_b
        def buddy(s):
            diff = abs(alpha - binom(s*p_b))
            hi = s + mpf(10**cur_pow10)
            hi_binom = binom(hi*p_b)
            if (abs(alpha-hi_binom)<diff):
                return buddy(hi)
            else:
                lo = s - mpf(10**cur_pow10)
                lo_binom = binom(lo*p_b)
                if (abs(alpha-lo_binom) < diff):
                    return buddy(lo)
                else:
                    return s
        for x in range(sigfigs):
            cur_val = buddy(cur_val)
            cur_pow10 -= 1
        if cur_val < 0:
            cur_val = 10**cur_pow10
        if cur_val > upper_bound:
            cur_val = upper_bound - 10**cur_pow10
        return round(cur_val,abs(int(cur_pow10)))
