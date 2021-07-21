from q_finder import *
from kardis_test import *
from data_binarizer import *
from s_prime import *
from mpmath import *

def hypothesis_test(data, value_list, alpha = 0.05, hypothesis = []):
    """
    Performs a multinomial SC test upon the dataset 'data'. 
    'value_list' is the list of values that may be observed within 'data', 
    'alpha' is the selected alpha level, and 'hypothesis' is a list representing 
    a probability distribution for each of the values in 'value_list'. 
    Returns the SC kardis, a boolean value for rejection, 
    as well as the lower plausibility bounds on s and p 
    and the closest plausible distribution if the given hypothesis is rejected.
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
        q = list(q_finder_main_slsqp(count_vector, hyp, p_lowerbound))
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
    Returns the same values as 'hypothesis_test'. 
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
    """
    Performs an exact binomial test (with tighter bounds than 'binary_hypothesis_test')
    upon the non-binary dataset 'data'. 'sigfigs' is the number of significant figures
    requested by the user, and all other inputs are the same as in 'binary_hypothesis_test'.
    Returns the closest plausible distribution, along with a boolean representing whether
    the original hypothesis was rejected. 
    """
    count_vector = [data.count(x) for x in value_list]
    s_prime, flipped = s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs)
    #the below print statements may need to be reworked
    if s_prime == 1:
        print("Proposed distribution not rejected at alpha = " + str(alpha) + ".")
        return (binary_hypothesis, False)
    elif flipped:
        q = [1-s_prime*binary_hypothesis[1], s_prime*binary_hypothesis[1]]
        print("Binomial tail exceeds 1 - "+ str(alpha) + ".")
    else:
        q = [s_prime*binary_hypothesis[0], 1-s_prime*binary_hypothesis[0]]
    q = [round(val, sigfigs-1) for val in q]
    print("Proposed distribution is rejected at alpha = " + str(alpha) + ". The closest plausible distribution is " + str(q) + ".")
    return (q, True)