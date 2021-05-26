import math

test_data = 20*[('male', 'white')] + [('female', 'white'), ('male', 'black'), ('female', 'asian')]

def multvarbiastest(data, numvars, alpha):
    var_dict = {}

    for i in range(numvars):
        var_dict[i] = []

    for i in range(len(data)):
        for j in range(numvars):
            if data[i][j] not in var_dict[j]:
                var_dict[j].append(data[i][j])
    '''
    for key in var_dict.keys():
        print(key, ' : ', var_dict[key])
    '''
    
    possible_combos = 1
    for key in var_dict.keys():
        possible_combos *= len(var_dict[key])

    norm_scriptx = possible_combos**len(data)
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))

    freq_dict = {}
    for entry in data:
        if entry in freq_dict:
            freq_dict[entry] += 1
        else:
            freq_dict[entry] = 1
    assumed_bias = max(freq_dict, key=freq_dict.get)
    #print('Assumed bias: ', assumed_bias)
    
    mg = 0
    for i in range(0, len(data)-data.count(assumed_bias)+1):
        mg += possible_combos * math.comb(len(data), i)
    
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