from Users.alessiaserafini.Desktop.amistadworkspace.q_finder_count_based import q_finder_main
from counts_kardis_test import *
from observation_count import *

"""
The observation and value_list should be in the format of a list of list
i.e for a 6-sided die:  
observation = 45*[[1]]+5*[[2]]+4*[[3]]+3*[[4]]+2*[[5]]+[[6]] 
value_list = [[1],[2],[3],[4],[5],[6]]
Output: [45, 5, 4, 3, 2, 1]
"""

def hypothesis_test(data, value_list, alpha, hypothesis = []):
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    count_vector = observation_count(data, value_list)
    kardis, reject, r, nu, h = kardis_test_main(count_vector, alpha, hypothesis)
    if reject:
        print("Proposed distribution rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        s_lowerbound = (alpha*nu)/(r*h)
        p_lowerbound = s_lowerbound*h
        print("Any plausible distribution must boost probability over the given distribution by " \
            + str(s_lowerbound) + ", and will therefore have a minimum probability of " + str(p_lowerbound) + ".")
        q = list(q_finder_main(count_vector, hypothesis, p_lowerbound))
        print("Closest plausible distribution: " + str(q))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        print("Proposed distribution not rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        return (kardis, reject)
    
