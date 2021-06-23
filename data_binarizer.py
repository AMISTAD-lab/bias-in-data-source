import math

def data_binarizer(nbdata, nbhypothesis, nbvaluelist, selectedvalue):
    s_index = nbvaluelist.index(selectedvalue)
    bvaluelist = [str(selectedvalue), 'not-'+str(selectedvalue)]
    bhypothesis = [nbhypothesis[s_index], math.fsum(nbhypothesis[:s_index]+nbhypothesis[s_index+1:])]
    bdata = [str(selectedvalue)]*nbdata.count(selectedvalue) 
    bdata += ['not-'+str(selectedvalue)]*(len(nbdata)-nbdata.count(selectedvalue))
    bcounts = [nbdata.count(selectedvalue), len(nbdata)-nbdata.count(selectedvalue)]
    return (bdata, bcounts, bhypothesis, bvaluelist)


def data_binarizer_multival(nbdata, nbhypothesis, nbvaluelist, selectedvaluelist):
    #probably a better way to get the name but. it's not coming to me rn
    biasval_name = '('+ str(selectedvaluelist[0])
    for val in selectedvaluelist[1:]:
        biasval_name += ", " + str(val)
    biasval_name+=')'
    bvaluelist = [biasval_name, 'not-'+biasval_name]
    num_biasval = sum([nbdata.count(x) for x in selectedvaluelist])
    bcounts = [num_biasval, len(nbdata)-num_biasval]
    bdata = [biasval_name]*num_biasval + ['not-'+biasval_name]*(len(nbdata)-num_biasval)
    bhypothesis = [math.fsum([nbhypothesis[nbvaluelist.index(x)] for x in selectedvaluelist])]
    bhypothesis += [1-bhypothesis[0]]
    return (bdata,bcounts,bhypothesis,bvaluelist)

def binarizer_sans_hyp(nbdata, selectedvaluelist):
    biasval_name = '('+ str(selectedvaluelist[0])
    for val in selectedvaluelist[1:]:
        biasval_name += ", " + str(val)
    biasval_name+=')'
    bvaluelist = [biasval_name, 'not-'+biasval_name]
    num_biasval = sum([nbdata.count(x) for x in selectedvaluelist])
    bcounts = [num_biasval, len(nbdata)-num_biasval]
    bdata = [biasval_name]*num_biasval + ['not-'+biasval_name]*(len(nbdata)-num_biasval)
    return (bdata,bcounts,bvaluelist)
    