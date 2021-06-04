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


def bin_grad(hyp,counts):
    return [counts[0]*(hyp[0]**(counts[0]-1))*(hyp[1]**counts[1]), \
        counts[1]*(hyp[0]**counts[0])*(hyp[1]**(counts[1]-1))]



