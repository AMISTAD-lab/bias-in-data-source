import math
from scipy.special import rel_entr

def kldiv(data, numvals=0, givendist = []):
    """calculates the kl of the observed distribution within data 
    to a given distribution. (if not given just uses uniform)
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
    if givendist == []:
        givendist = [1/len(proplist)]*len(proplist)
    elif not (len(givendist) == numvals):
        raise Exception("given distribution mismatch with values")
    kl = 0
    for index in range(len(proplist)):
        kl += proplist[index]*math.log(proplist[index]/givendist[index])
    return kl

def simplekldiv(p,q):
    """wherein both p and q are distributions"""
    kl = 0
    for index in range(len(p)):
        kl+= p[index]*math.log(p[index]/q[index])
    return kl

def trinarybruteforcemg(size, min_kl, givendist = []):
    """i want to extend this to arbitrary nb situations but trinary first
    bc this might need recursion and i don't want recursion"""
    mg = 0
    if givendist == []:
        givendist = [1/3]*3 #uniform dist
    for q in range(size+1):   
        for r in range(size-q+1):
            s=size-q-r
            #print(q,r,s)
            #plist = [q/size, r/size, s/size]
            #kl = sum(rel_entr(plist,givendist))
            fakedata = [0]*q + [1]*r + [2]*s
            kl = kldiv(fakedata, numvals=3, givendist=givendist)
            #print(kl)
            if kl >= min_kl:
                #print(q,r,s)
                combs = math.comb(size,q)*math.comb(size-q,r)*math.comb(size-q-r,s)
                mg += combs
    #this is absolutely going to need recursion i am suffering
    return mg

def trinarymgcalculator(data, givendist = []):
    min_kl = kldiv(data, numvals=3, givendist=givendist)
    #print("minkl",min_kl)
    mg = trinarybruteforcemg(len(data), min_kl, givendist=givendist)
    return mg


            



