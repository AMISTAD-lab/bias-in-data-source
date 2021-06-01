import math

# Test inputs for a hiring process
test_data = 20*[('male', 'white')] + [('female', 'white'), ('male', 'black'), ('female', 'asian')]
test_numvars = 2
test_varpossiblevals = [2, 3]
test_assumed_bias = ('male', 'white')
test_alpha = 0.05

# Test inputs for rolling a die 20 times
test_data2 = 15*[(1)] + [(2), (3), (4), (5), (6)]
test_numvars2 = 1
test_varpossiblevals2 = [6]
test_assumed_bias2 = (1)
test_alpha2 = 0.05

def multvarbiastest(data, numvars, varpossiblevals, assumed_bias, alpha):
    """
    Notes:
    1. data should be a list of tuples, not lists
    2. varpossiblevalues is a list of possible values for each variable
       in the order in which they appear in the tuples of data
    """

    # Computes all possible combinations for a given data point
    possible_combos = 1
    for i in range(numvars):
        possible_combos *= varpossiblevals[i]

    norm_scriptx = possible_combos**len(data)
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))

    # Assumed bias is the most common data entry
    """
    freq_dict = {}
    for entry in data:
        if entry in freq_dict:
            freq_dict[entry] += 1
        else:
            freq_dict[entry] = 1
    assumed_bias = max(freq_dict, key=freq_dict.get)
    """
    bias = assumed_bias
    
    mg = 0
    for i in range(0, len(data)-data.count(bias)+1):
        mg += possible_combos * math.comb(len(data), i) * (possible_combos-1)**i
    
    nu = norm_scriptx/mg
    kardis = r*u/nu
    s_lowerbound = alpha*nu/(r*u)
    if kardis < alpha:
        reject = True
        print("Uniform distribution rejected at alpha = "+ str(alpha) +". Kardis = " + str(kardis)+ ". s >", s_lowerbound)
    else:
        reject = False
        print("Uniform distribution not rejected at alpha = "+ str(alpha) +". Kardis = " + str(kardis)+ ". s >", s_lowerbound)
    return (kardis, s_lowerbound, reject)