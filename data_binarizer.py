import math

def data_binarizer(nbdata, nbhypothesis, nbvaluelist, selectedvalue):
    s_index = nbvaluelist.index(selectedvalue)
    bvaluelist = [str(selectedvalue), 'not-'+str(selectedvalue)]
    bhypothesis = [nbhypothesis[s_index], math.fsum(nbhypothesis[:s_index]+nbhypothesis[s_index+1:])]
    bdata = [str(selectedvalue)]*nbdata.count(selectedvalue) 
    bdata += ['not-'+str(selectedvalue)]*(len(nbdata)-nbdata.count(selectedvalue))
    bcounts = [nbdata.count(selectedvalue), len(nbdata)-nbdata.count(selectedvalue)]
    return (bdata, bcounts, bhypothesis, bvaluelist)