from scipy import stats
import math
from data_binarizer import *
from mpmath import *

def s_prime_finder(data, hypothesis, value_list, bias_value_list, alpha, sigfigs):
    """ Returns the lower bound on the scalar that your original hypothesis value 
        for bias_value must be multiplied by in order to return a plausible explanation 
        at your given alpha. for instance, 
        s_prime_finder(10*[0]+5*[1], [0.2,0.8], [0,1], 0, 0.05, 4) 
        tells you that 0.2 is too low, and must be multiplied by 2.1957 in order
        to be a reasonable probability for the value 0.
    """
    bdata, bcounts, bhyp, bvaluelist = data_binarizer_multival(data,hypothesis,value_list,bias_value_list)
    assert bhyp[0] + bhyp[1] == 1
    n = len(data)
    i = bcounts[0]
    def stirl_approx(k):
        return (mpf(n/math.e)**n)/(sqrt(2*math.pi*k*(n-k)) * mpf(k/math.e)**k * mpf((n-k)/math.e)**(n-k))
    prob_more_extreme = sum([stirl_approx(k)*(mpf(bhyp[0])**k)*(mpf(bhyp[1])**(n-k)) for k in range(i,n)])
    prob_more_extreme += mpf(bhyp[0])**n
    #print("prob gotten")
    if prob_more_extreme >= alpha:
        return 1 #your explanation is fine.
    else:
        s_lowerbound = alpha/prob_more_extreme
        target = s_lowerbound* (mpf(bhyp[0])**i * mpf(bhyp[1])**(n-i))
        cur_val = mpf(s_lowerbound)**(1/mpf(bcounts[0]))
        #print(cur_val)
        starter_pow10 = math.log10(cur_val)//1
        cur_pow10 = starter_pow10
        #print("time for recursion")
        while cur_pow10 > starter_pow10-sigfigs:
            #print(cur_pow10)
            cur_val = buddy(cur_pow10,cur_val,n,i,bhyp[0],target)
            cur_pow10 -= 1
        #print(cur_val)
        return round(cur_val,abs(int(cur_pow10)))

def s_prime_finder2(nbdata, bhyp, nbvalue_list, bias_value_list, alpha, sigfigs):
    """for giving a binary hypothesis instead of an nb one."""
    bdata,bcounts,bvaluelist = binarizer_sans_hyp(nbdata, bias_value_list)
    n = len(nbdata)
    i = bcounts[0]
    def stirl_approx(k):
        return (mpf(n/math.e)**n)/(sqrt(2*math.pi*k*(n-k)) * mpf(k/math.e)**k * mpf((n-k)/math.e)**(n-k))
    prob_more_extreme = sum([stirl_approx(k)*(mpf(bhyp[0])**k)*(mpf(bhyp[1])**(n-k)) for k in range(i,n)])
    prob_more_extreme += mpf(bhyp[0])**n
    #print("prob gotten")
    if prob_more_extreme >= alpha:
        return 1 #your explanation is fine.
    else:
        s_lowerbound = alpha/prob_more_extreme
        target = s_lowerbound* (mpf(bhyp[0])**i * mpf(bhyp[1])**(n-i))
        cur_val = mpf(s_lowerbound)**(1/mpf(bcounts[0]))
        #print(cur_val)
        starter_pow10 = math.log10(cur_val)//1
        cur_pow10 = starter_pow10
        #print("time for recursion")
        while cur_pow10 > starter_pow10-sigfigs:
            #print(cur_pow10)
            cur_val = buddy(cur_pow10,cur_val,n,i,bhyp[0],target)
            cur_pow10 -= 1
        #print(cur_val)
        return round(cur_val,abs(int(cur_pow10)))

def s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs):
    bcounts = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(count_vector)
    i = bcounts[0]
    def stirl_approx(k):
        return (mpf(n/math.e)**n)/(sqrt(2*math.pi*k*(n-k)) * mpf(k/math.e)**k * mpf((n-k)/math.e)**(n-k))
    prob_more_extreme = sum([stirl_approx(k)*(mpf(binary_hypothesis[0])**k) \
        *(mpf(binary_hypothesis[1])**(n-k)) for k in range(i,n)])
    prob_more_extreme += mpf(binary_hypothesis[0])**n
    if prob_more_extreme >= alpha:
        return 1 #your explanation is fine.
    else:
        s_lowerbound = alpha/prob_more_extreme
        target = s_lowerbound* (mpf(binary_hypothesis[0])**i * mpf(binary_hypothesis[1])**(n-i))
        cur_val = mpf(s_lowerbound)**(1/mpf(bcounts[0]))
        starter_pow10 = math.log10(cur_val)//1
        cur_pow10 = starter_pow10
        while cur_pow10 > starter_pow10-sigfigs:
            #print(cur_pow10)
            cur_val = buddy(cur_pow10,cur_val,n,i,binary_hypothesis[0],target)
            cur_pow10 -= 1
        #print(cur_val)
        return round(cur_val,abs(int(cur_pow10)))

def buddy(cur_pow10,cur_val,n,i,pb,target):
    """recursively optimizes cur_val for a single power of 10"""
    cur_prob = (cur_val*mpf(pb))**i * (1-cur_val*mpf(pb))**(n-i)
    diff = abs(target - cur_prob)
    lo = cur_val - mpf(10**cur_pow10)
    lo_prob = (lo*mpf(pb))**i * (1-lo*mpf(pb))**(n-i)
    if abs(target-lo_prob) < diff:
        return buddy(cur_pow10,lo,n,i,pb,target)
    else:
        hi = cur_val + 10**cur_pow10
        hi_prob = (hi*mpf(pb))**i * (1-hi*mpf(pb))**(n-i)
        if abs(target-hi_prob) < diff:
            return buddy(cur_pow10,hi,n,i,pb,target)
        else:
            return cur_val