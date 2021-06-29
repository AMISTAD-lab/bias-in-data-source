from q_finder_count_based import *
from counts_kardis_test import *
from data_binarizer import *

def hypothesis_test(data, value_list, alpha = 0.05, hypothesis = []):
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    count_vector = [data.count(x) for x in value_list]
    kardis, reject, r, nu, h = kardis_test_main(count_vector, alpha, hyp)
    if reject:
        print("Proposed distribution rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        s_lowerbound = (alpha*nu)/(r*h)
        p_lowerbound = s_lowerbound*h
        print("Any plausible distribution must boost probability over the given distribution by " \
            + str(s_lowerbound) + ", and will therefore have a minimum probability of " + str(p_lowerbound) + ".")
        q = list(q_finder_main(count_vector, hyp, p_lowerbound))
        print("Closest plausible distribution: " + str(q))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        print("Proposed distribution not rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        return (kardis, reject)

def binary_hypothesis_test(data, value_list, selected_value_list, alpha=0.05, binary_hypothesis=[0.5,0.5]):
    count_vector = [data.count(x) for x in value_list]
    binary_count_vector = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(binary_count_vector)
    k = binary_count_vector[0]
    kardis, reject, r, nu, h = kardis_test_main(binary_count_vector, alpha, binary_hypothesis)
    if reject:
        print("Proposed distribution rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        s_lowerbound = alpha * ((n+1) * math.comb(n,k) * binary_hypothesis[0]**k * binary_hypothesis[1]**(n-k))**(-1)
        p_lowerbound = s_lowerbound*h
        print("Any plausible distribution must boost probability over the given distribution by " \
            + str(s_lowerbound) + ", and will therefore have a minimum probability of " + str(p_lowerbound) + ".")
        q = list(binary_q_finder_main(binary_count_vector, binary_hypothesis, alpha))
        print("Closest plausible distribution: " + str(q))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        print("Proposed distribution not rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        return (kardis, reject)