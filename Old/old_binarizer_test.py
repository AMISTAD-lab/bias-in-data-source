from q_finder_count_based import *
from counts_kardis_test import *
from s_prime import *
from data_binarizer import *

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