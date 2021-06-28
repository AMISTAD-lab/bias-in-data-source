from q_finder_count_based import *
from counts_kardis_test import *
from s_prime import *

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

def binarizer_test(data, value_list, selected_value_list, alpha = 0.05, binary_hypothesis = [0.5,0.5], sigfigs = 3):
    count_vector = [data.count(x) for x in value_list]
    s_prime, flipped = s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs)
    #the below print statements may need to be reworked
    if s_prime == 1:
        print("Proposed distribution not rejected at alpha = " + str(alpha) + ".")
    elif flipped:
        print("Binomial tail exceeds 1 - "+ str(alpha)+". Proposed probability of non-selected values must be multiplied by "\
            +str(s_prime)+" to become a valid explanation")
    else:
        print("Proposed probability of selected values must be multiplied by "\
            +str(s_prime)+" to become a valid explanation.")
    return s_prime, flipped
