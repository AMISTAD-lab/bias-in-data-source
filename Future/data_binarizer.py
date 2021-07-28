import math

def data_binarizer_main(nbcounts, nbvalue_list, selected_value_list):
    """
    Binarizes the provided non-binary count vector 'nbcounts'.
    'nbvalue_list' is the list of values found in the data, 
    where nbcounts[i] contains the number of times nbvalue_list[i]
    occurs in a dataset. Returns the binary count vector containing
    the quantity of values found in 'selected_value_list' 
    and the quantity of not-selected values.
    """
    num_biasval = sum([nbcounts[nbvalue_list.index(x)] for x in selected_value_list])
    bcounts = [num_biasval, sum(nbcounts)-num_biasval]
    return bcounts
    
def data_binarizer(nbdata, nbhypothesis, nbvaluelist, selectedvalue):
    """
    A deprecated data binarization function.
    """
    s_index = nbvaluelist.index(selectedvalue)
    bvaluelist = [str(selectedvalue), 'not-'+str(selectedvalue)]
    bhypothesis = [nbhypothesis[s_index], math.fsum(nbhypothesis[:s_index]+nbhypothesis[s_index+1:])]
    bdata = [str(selectedvalue)]*nbdata.count(selectedvalue) 
    bdata += ['not-'+str(selectedvalue)]*(len(nbdata)-nbdata.count(selectedvalue))
    bcounts = [nbdata.count(selectedvalue), len(nbdata)-nbdata.count(selectedvalue)]
    return (bdata, bcounts, bhypothesis, bvaluelist)

def data_binarizer_multival(nbdata, nbhypothesis, nbvaluelist, selectedvaluelist):
    """
    A deprecated data binarization function. 
    """
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
    """
    A deprecated data binarization function. 
    """
    biasval_name = '('+ str(selectedvaluelist[0])
    for val in selectedvaluelist[1:]:
        biasval_name += ", " + str(val)
    biasval_name+=')'
    bvaluelist = [biasval_name, 'not-'+biasval_name]
    num_biasval = sum([nbdata.count(x) for x in selectedvaluelist])
    bcounts = [num_biasval, len(nbdata)-num_biasval]
    bdata = [biasval_name]*num_biasval + ['not-'+biasval_name]*(len(nbdata)-num_biasval)
    return (bdata,bcounts,bvaluelist)
