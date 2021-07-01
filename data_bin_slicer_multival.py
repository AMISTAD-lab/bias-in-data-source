def data_bin_slicer_multival(nbdata, nbhypothesis, nbvaluelist, selectedvaluelist):
    """
    	Similar to data_binarizer,except the selected values are not combined to 1 variable
	ex. nbvaluelist = [1,2,3,4]
	selectedvaluelist = [1,2]
	nvaluelist = [[1],[2],["not-([1,2])"]]
	"""
    nvaluelist = []
   
    nhypothesis = []
    ndata = []
    biasval_name = '('+ str(selectedvaluelist[0])
    for val in selectedvaluelist[1:]:
        biasval_name += ", " + str(val)
    biasval_name+=')'
    
    for val in range(len(selectedvaluelist)):
        nvaluelist.append(selectedvaluelist[val])
    
    nvaluelist.append(['not-'+biasval_name])
    num_biasval = [nbdata.count(x) for x in selectedvaluelist]
    ncounts = num_biasval
    
    ncounts.append(int(len(nbdata))-int(sum(num_biasval)))
    for x in range(len(selectedvaluelist)):
        ndata += [selectedvaluelist[x]]*num_biasval[x]
 
    ndata += [['not-'+biasval_name]]*ncounts[-1]
            
    for y in range(len(selectedvaluelist)):
        nhypothesis.append(nbhypothesis[nbvaluelist.index(selectedvaluelist[y])])
        
    sum_nhyp = 0
   
    for s in range(len(nhypothesis)):
        sum_nhyp += nhypothesis[s]
    nhypothesis.append(1-sum_nhyp)
    return (ndata,ncounts,nhypothesis,nvaluelist)
