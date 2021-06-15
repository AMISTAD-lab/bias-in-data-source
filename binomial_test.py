from scipy import stats

"""
Example tested:

data = 45*[1] + 5*[2] + 4*[3] + 3*[4] + 2*[5] + [6]
hyp = [1/6]*6
value_list = [1,2,3,4,5,6]
bias = 1

"""

def binomial_test(data,hypothesis,value_list,bias_value):
    bdata, bcounts, bhypothesis, bvaluelist = data_binarizer(data,hypothesis,value_list,bias_value)
    num_events = len(data)
    idx = value_list.index(int(bvaluelist[0])) #gets index of the assumed bias value
    value_prob = hypothesis[idx] 			   #gives the original hypothesis for the assumed bias value
    p_value = stats.binom_test(bcounts, n=num_events, p=value_prob, alternative='greater')
    return p_value