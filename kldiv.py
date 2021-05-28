import math

def kldiv_fromuniform(data, numvals=0):
    """calculates the kl of the observed distribution within data to a uniform dist.
    set numvals to number of possible values w/in data if not all are observed"""
    valuedict = {}
    for item in data:
        if item not in valuedict.keys():
            valuedict[item] = 1
        else:
            valuedict[item] +=1
    if numvals == 0:
        numvals = len(valuedict.keys())
    proplist = [valuedict[x]/len(data) for x in valuedict.keys()]
    #for values not observed in data
    #cannot actually be 0, so we approximate
    proplist += [0.000000001]*(numvals-len(valuedict.keys())) 
    uniproplist = [1/len(proplist)]*len(proplist)
    kl = 0
    for index in range(len(proplist)):
        kl += proplist[index]*math.log2(proplist[index]/uniproplist[index])
    return kl

def simplekldiv(p,q):
    """wherein both p and q are distributions"""
    kl = 0
    for index in range(len(p)):
        kl+= p[index]*math.log2(p[index]/q[index])
    return kl
