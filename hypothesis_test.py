from q_finder import *
from kardis_test import *
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

def hypothesis_test_silent(data, value_list, alpha = 0.05, hypothesis = []):
    """
    Same as `hypothesis_test`, but without print statements
    """
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    count_vector = [data.count(x) for x in value_list]
    kardis, reject, r, nu, h = kardis_test_main(count_vector, alpha, hyp)
    if reject:
        s_lowerbound = (alpha*nu)/(r*h)
        p_lowerbound = s_lowerbound*h
        q = list(q_finder_main_slsqp(count_vector, hyp, p_lowerbound))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        return (kardis, reject)