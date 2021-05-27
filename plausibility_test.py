import math
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


def givennumberplausibility(data,alpha,numberofvalues,hyp,hyp_var):
    kardis, s_lowerbound, p_lowerbound, reject = givennumberbiastest(data, numberofvalues, alpha)
    norm_scriptx = (numberofvalues**(len(data)))
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
