import math

sampledata = [0,0,0,1,1,0,0,1,0,0,1]

def binarybiastest(data,alpha):
    """assumes data is a binary list."""
    norm_scriptx = (2**(len(data)))
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))
    if data.count(0) > data.count(1):
        assumed_bias = 0
    else:
        assumed_bias = 1
    mg = 0
    for i in range(0,len(data)-data.count(assumed_bias)+1):
        mg += math.comb(len(data), i)
    mg = mg*2
    nu = norm_scriptx/mg
    #print(nu)
    kardis = r*u/nu
    s_lowerbound = alpha*nu/(r*u)
    if kardis < alpha:
        reject = True
        print("Uniform distribution rejected at alpha = "+ str(alpha) +". Kardis = " + str(kardis)+ ". s >", s_lowerbound)
    else:
        reject = False
        print("Uniform distribution not rejected at alpha = "+ str(alpha) +". Kardis = " + str(kardis)+ ". s >", s_lowerbound)
    
    return (kardis, s_lowerbound, reject)
    
def alltypebiastest(data,alpha):
    """takes in a list of values and does bias testing on them.
        Assumes that the list contains all potential values of data."""
    str_to_int_dict = {}
    int_list_representation = []
    #the below code takes in your whatever-type list and represents it as a list of ints
    #if you have 3 different values in your list, the output will contain the values 0, 1, 2
    #the dictionary tells you which output value corresponds to the input value
    for element in data:
        if element not in str_to_int_dict.keys():
            str_to_int_dict[element] = len(str_to_int_dict.keys())
        int_list_representation += [str_to_int_dict[element]] 
    #calculates norm_scriptx based on number of unique values in data, raised to its length
    norm_scriptx = (len(str_to_int_dict.keys())**(len(data)))
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))
    #elecount contains the number of times a value shows up in the list,
    #at the index of its corresponding integer.
    #so if m corresponds to 0, and m shows up 5 times,
    #elecount[0] = 5
    elecount = [int_list_representation.count(x) for x in str_to_int_dict.values()]
    #assumed_bias is the integer representation of the value
    #that shows up the most in the data
    assumed_bias = elecount.index(max(elecount))
    #this variable is the original value that we're assuming the data is biased towards
    assumed_bias_value = list(str_to_int_dict.keys())[list(str_to_int_dict.values()).index(assumed_bias)]
    mg = 0
    #this sets mg to the number of sequences that are at least as biased as our sequence
    for i in range(0,len(data)-int_list_representation.count(assumed_bias)+1):
        mg += math.comb(len(data),i)*((len(str_to_int_dict.keys())-1)**i)
    mg = mg*len(str_to_int_dict.keys())
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

def givennumberbiastest(data, numberofvalues, alpha):
    """takes in a list of values and does bias testing on them.
        you input how many potential values your data may have (like 2 if binary)
        please make sure you don't put in less values than do exist"""
    str_to_int_dict = {}
    int_list_representation = []
    #the below code takes in your whatever-type list and represents it as a list of ints
    #if you have 3 different values in your list, the output will contain the values 0, 1, 2
    #the dictionary tells you which output value corresponds to the input value
    for element in data:
        if element not in str_to_int_dict.keys():
            str_to_int_dict[element] = len(str_to_int_dict.keys())
        int_list_representation += [str_to_int_dict[element]] 
    #calculates norm_scriptx based on the given # of values in data, raised to its length
    norm_scriptx = (numberofvalues**(len(data)))
    u = 1/norm_scriptx
    r = norm_scriptx*(1+math.log(norm_scriptx))
    #elecount contains the number of times a value shows up in the list,
    #at the index of its corresponding integer.
    #so if m corresponds to 0, and m shows up 5 times,
    #elecount[0] = 5
    elecount = [int_list_representation.count(x) for x in str_to_int_dict.values()]
    #assumed_bias is the integer representation of the value
    #that shows up the most in the data
    assumed_bias = elecount.index(max(elecount))
    #this variable is the original value that we're assuming the data is biased towards
    assumed_bias_value = list(str_to_int_dict.keys())[list(str_to_int_dict.values()).index(assumed_bias)]
    mg = 0
    #this sets mg to the number of sequences that are at least as biased as our sequence
    for i in range(0,len(data)-int_list_representation.count(assumed_bias)+1):
        mg += math.comb(len(data),i)*((numberofvalues-1)**i)
    mg = mg*numberofvalues
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
