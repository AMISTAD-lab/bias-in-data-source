from q_finder_count_based import *
from counts_kardis_test import *
from data_binarizer import *
from s_prime import *
from mpmath import *

def hypothesis_test(data, value_list, alpha = 0.05, hypothesis = []):
    """
    Performs a multinomial SC test upon the dataset 'data'. 
    'value_list' is the list of values that may be observed within 'data', 
    'alpha' is the selected alpha level, and 'hypothesis' is a list of probability values
    in the range (0,1] where hypothesis[i] corresponds to the probability of observing
    the value in value_list[i]. Returns the SC kardis, a boolean value for rejection, 
    as well as the lower plausibility bounds on s and p and the closest plausible distribution q
    if the given hypothesis is rejected.
    """
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
    """
    Performs a binomial SC test upon the non-binary dataset 'data'. 
    'selected_value_list' is a list of the values the user believes the data is biased towards. 
    The data is transformed into a binary dataset by use of that list, 
    as a collection of selected values and not-selected values. 
    'binary_hypothesis' contains the probability of producing a selected value,
    and the probability of producing a not-selected value.
    'value_list' and 'alpha' are the same as in 'hypothesis_test'.
    
    """
    count_vector = [data.count(x) for x in value_list]
    binary_count_vector = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(binary_count_vector)
    k = binary_count_vector[0]
    kardis, reject, r, nu, h = binary_kardis_test_main(binary_count_vector, alpha, binary_hypothesis)
    if reject:
        print("Proposed distribution rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        s_lowerbound = alpha * ((n+1) * math.comb(n,k)* mpf(binary_hypothesis[0])**k * mpf(binary_hypothesis[1])**(n-k))**(-1)
        p_lowerbound = s_lowerbound*h
        print("Any plausible distribution must boost probability over the given distribution by " \
            + str(s_lowerbound) + ", and will therefore have a minimum probability of " + str(p_lowerbound) + ".")
        q = list(binary_q_finder_main(binary_count_vector, binary_hypothesis, alpha))
        print("Closest plausible distribution: " + str(q))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        print("Proposed distribution not rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        return (kardis, reject)
  
def exact_binomial_test(data, value_list, selected_value_list, alpha = 0.05, binary_hypothesis = [0.5,0.5], sigfigs = 4):
    count_vector = [data.count(x) for x in value_list]
    s_prime, flipped = s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs)
    #the below print statements may need to be reworked
    if s_prime == 1:
        print("Proposed distribution not rejected at alpha = " + str(alpha) + ".")
    elif flipped:
        print("Binomial tail exceeds 1 - "+ str(alpha)+ \
            ". Proposed probability of non-selected values must be multiplied by "\
            +str(s_prime)+" to become a valid explanation")
    else:
        print("Proposed probability of selected values must be multiplied by "\
            +str(s_prime)+" to become a valid explanation.")
    return s_prime, flipped
