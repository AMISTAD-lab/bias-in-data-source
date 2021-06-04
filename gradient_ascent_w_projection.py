import math
import decimal
#sample data and hypothesis for short times:
"""
data = 35*[0]+2*[1]
hyp = [0.7,0.3]
sux = 1.3326866953238921*10**(-6)
"""
#yes step size is huge but gradient is very very small
def binary_grad_ascent(data, hypothesis, sux, step_size = 100000000):
    """does gradient ascent on probability, but this time to maintain the sum-to-one constraint
    I project each stepped value back onto the line y=1-x instead of making a penalty function"""
    counts = [decimal.Decimal(data.count(0)), decimal.Decimal(data.count(1))]
    hyp = [decimal.Decimal(hypothesis[0]),decimal.Decimal(hypothesis[1])]
    step_size=decimal.Decimal(step_size)
    while (hyp[0]**counts[0])*(hyp[1]**counts[1])<sux:
        gradient = bin_grad(hyp,counts)
        temphyp = [hyp[0]+gradient[0], hyp[1]+gradient[1]]
        hypval = (temphyp[0]-temphyp[1]+1)/2
        hyp = [hypval,1-hypval]
        print(hyp)
    return hyp

def grad_ascent(data,value_list,hypothesis,sux,step_size = 100000000):
    counts = [decimal.Decimal(data.count(val)) for val in value_list]
    hyp = [decimal.Decimal(hypval) for hypval in hypothesis]
    step_size=decimal.Decimal(step_size)
    while nb_probability_func(hyp,counts) < sux:
        gradient = nb_probability_grad(hyp,counts)
        temphyp = [hyp[i]+gradient[i] for i in range(len(hyp))]
        norm_vector_scalar = (1-sum(temphyp))/len(temphyp)
        hyp = [temphyp[i] + norm_vector_scalar for i in range(len(temphyp))]
    return hyp

def nb_probability_func(hyp,counts):
    func = 1
    for i in range(len(hyp)):
        func = func*(hyp[i]**counts[i])
    return func

def nb_probability_grad(hyp,counts):
    grad = []
    for i in range(len(hyp)):
        hypslice = hyp[:i]+hyp[i+1:]
        countsslice = counts[:i]+counts[i+1:]
        grad += [counts[i]*(hyp[i]**(counts[i]-1))*nb_probability_func(hypslice,countsslice)]
    return grad


def bin_grad(hyp,counts):
    return [counts[0]*(hyp[0]**(counts[0]-1))*(hyp[1]**counts[1]), \
        counts[1]*(hyp[0]**counts[0])*(hyp[1]**(counts[1]-1))]



