from q_finder_count_based import *
from counts_kardis_test import *
from s_prime import *

def binary_hypothesis_test(data, alpha = 0.05, hypothesis = []):
    if hypothesis == []:
        hyp = [0.5, 0.5]
    else:
        hyp = hypothesis
    count_vector = [data.count(1), data.count(0)]
    n = sum(count_vector)
    k = count_vector[0]
    kardis, reject, r, nu, h = kardis_test_main(count_vector, alpha, hyp)
    if reject:
        print("Proposed distribution rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        s_lowerbound = alpha * ((n+1) * math.comb(n,k) * hyp[0]**k * hyp[1]**(n-k))**(-1)
        p_lowerbound = s_lowerbound*h
        print("Any plausible distribution must boost probability over the given distribution by " \
            + str(s_lowerbound) + ", and will therefore have a minimum probability of " + str(p_lowerbound) + ".")
        q = list(binary_q_finder_main(count_vector, hyp, alpha))
        print("Closest plausible distribution: " + str(q))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        print("Proposed distribution not rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        return (kardis, reject)