import math
import numpy as np
import copy
import decimal

"""VERY broken at the moment."""

def binary_grad_ascent(data, hypothesis, sux, penalty_scalar = 1000, step_size = 0.01):
    """for data of 0s and 1s."""
    counts = [decimal.Decimal(data.count(0)), decimal.Decimal(data.count(1))]
    hyp = [decimal.Decimal(hypothesis[0]),decimal.Decimal(hypothesis[1])]
    step_size=decimal.Decimal(step_size)
    #p(x) function = (hyp[0]**counts[0])*(hyp[1]**counts[1])
    #penalty function = -penalty_scalar*(hyp[0]+hyp[1]-1)**2
    #these are added together to get f(x) which we maximize
    scalar = decimal.Decimal(penalty_scalar)
    while (hyp[0]**counts[0])*(hyp[1]**counts[1])<sux:
        gradient = bin_grad(hyp,counts,scalar)
        #print(gradient)
        hyp[0] = hyp[0] + step_size*gradient[0]
        hyp[1] = hyp[1] + step_size*gradient[1]
        print(hyp)
    return hyp

def func(hyp,counts,scalar):
    return (hyp[0]**counts[0])*(hyp[1]**counts[1])-scalar*(hyp[0]+hyp[1]-1)**2

def bin_grad(hyp,counts, scalar):
    return [counts[0]*(hyp[0]**(counts[0]-1))*(hyp[1]**counts[1]) - 2*scalar*(hyp[0]+hyp[1]-1),\
        (hyp[0]**counts[0])*counts[1]*(hyp[1]**(counts[1]-1))-2*scalar*(hyp[0]+hyp[1]-1)]
    