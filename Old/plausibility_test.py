import math
from biastest import *

s_lowerbound = 1764
data = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
alpha = 0.05
hyp = 9/10
def binaryplausibility(data,alpha,hyp,hyp_var):
    #hyp_var is the variable we're testing the hypothesis for
    kardis, s_lowerbound, reject = binarybiastest(data,alpha)
    norm_scriptx = (2**(len(data)))
    l = len(data)
    mechanism = hyp**data.count(hyp_var)*(1-hyp)**(l-data.count(hyp_var))
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))
    sux = s_lowerbound*u
    q = sux**(1/l)
    print(q)

    
    if hyp < q:
        plausible = False
    else:
        plausible = True
    
    return plausible


def givennumberplausibility(data,alpha,numberofvalues,hyp_var,var_prob_dict):
    """var_prob is a dictionary of the given probability for each variable to occur."""
    kardis, s_lowerbound, p_lowerbound, reject = givennumberbiastest(data, numberofvalues, alpha)
    mechanism = 1     
    norm_scriptx = (numberofvalues**(len(data)))
    l = len(data)
    if hyp_var in var_prob_dict.keys(): 
        hyp = var_prob_dict[hyp_var]
    else:
        hyp = data.count(hyp_var)/len(data)
    
    for x in var_prob_dict.keys():
        mechanism = mechanism*(var_prob_dict[x]**(data.count(x)))
    #mechanism = hyp**data.count(hyp_var)*(1-hyp)**(l-data.count(hyp_var))
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))
    sux = s_lowerbound*u
    q = sux**(1/l)
    #print(q)

    
    if hyp < q:
        plausible = False
    else:
        plausible = True
    
    return plausible
