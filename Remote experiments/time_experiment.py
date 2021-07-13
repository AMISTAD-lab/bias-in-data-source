import pandas as pd
import math
from itertools import *
import timeit
from mpmath import *
import numpy as np
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, Bounds
from scipy.special import rel_entr
import warnings
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager

def excel_to_list(filepath):
    data = pd.read_excel(filepath) 
    df = pd.DataFrame(data)
    array = df.to_numpy()
    data_list = array.tolist()
    return data_list
def excel_to_list(filepath,columns):
    data = pd.read_excel(filepath) 
    df = pd.DataFrame(data,columns= columns)
    array = df.to_numpy()
    data_list = array.tolist()
    return data_list

filepath = r"compass_100k.xlsx"
columns= ['sex','race']
observation = excel_to_list(filepath,columns) 
#print(len(observation))
#print(listg)
#print(type(observation))
if len(observation) == 100000 and type(observation) == list:
    print('Data processed successfully.')
else:
    print('Data processing failed.')

#----------------------------------------------------------------------------------------------

#from q_finder_count_based import *
#from counts_kardis_test import *
#from data_binarizer import *
#from s_prime import *

def hypothesis_test(data, value_list, alpha = 0.05, hypothesis = []):
    """
    Performs a multinomial SC test upon the dataset 'data'. 
    'value_list' is the list of values that may be observed within 'data', 
    'alpha' is the selected alpha level, and 'hypothesis' is a list representing 
    a probability distribution for each of the values in 'value_list'. 
    Returns the SC kardis, a boolean value for rejection, 
    as well as the lower plausibility bounds on s and p 
    and the closest plausible distribution if the given hypothesis is rejected.
    """
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    count_vector = [data.count(x) for x in value_list]
    kardis, reject, r, nu, h = kardis_test_main(count_vector, alpha, hyp)
    if reject:
        #print("Proposed distribution rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        s_lowerbound = (alpha*nu)/(r*h)
        p_lowerbound = s_lowerbound*h
        #print("Any plausible distribution must boost probability over the given distribution by " \
            #+ str(s_lowerbound) + ", and will therefore have a minimum probability of " + str(p_lowerbound) + ".")
        q = list(q_finder_main(count_vector, hyp, p_lowerbound))
        #print("Closest plausible distribution: " + str(q))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        #print("Proposed distribution not rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        return (kardis, reject)

def binary_hypothesis_test(data, value_list, selected_value_list, alpha=0.05, binary_hypothesis=[0.5,0.5]):
    """
    Performs a binomial SC test upon the non-binary dataset 'data'. 
    'selected_value_list' is a list of the values the user believes the data is biased towards. 
    The data is transformed into a binary dataset by use of that list, 
    as a collection of selected values and not-selected values. 
    'binary_hypothesis' contains the probability of producing a selected value,
    and the probability of producing a not-selected value.
    'value_list' and 'alpha' are the same as in 'hypothesis_test'.
    Returns the same values as 'hypothesis_test'. 
    """
    count_vector = [data.count(x) for x in value_list]
    binary_count_vector = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(binary_count_vector)
    k = binary_count_vector[0]
    kardis, reject, r, nu, h = binary_kardis_test_main(binary_count_vector, alpha, binary_hypothesis)
    if reject:
        print("Proposed distribution rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        s_lowerbound = alpha * ((n+1) * math.comb(n,k)* mpf(binary_hypothesis[0])**k * mpf(binary_hypothesis[1])**(n-k))**(-1)
        p_lowerbound = s_lowerbound*h
        print("Any plausible distribution must boost probability over the given distribution by " \
            + str(s_lowerbound) + ", and will therefore have a minimum probability of " + str(p_lowerbound) + ".")
        q = list(binary_q_finder_main(binary_count_vector, binary_hypothesis, alpha))
        print("Closest plausible distribution: " + str(q))
        return (kardis, reject, s_lowerbound, p_lowerbound, q)
    else:
        print("Proposed distribution not rejected at alpha = " + str(alpha)  + ". Kardis = " + str(kardis) + ".")
        return (kardis, reject)


def exact_binomial_test(data, value_list, selected_value_list, alpha = 0.05, binary_hypothesis = [0.5,0.5], sigfigs = 4):
    """
    Performs an exact binomial test (with tighter bounds than 'binary_hypothesis_test')
    upon the non-binary dataset 'data'. 'sigfigs' is the number of significant figures
    requested by the user, and all other inputs are the same as in 'binary_hypothesis_test'.
    Returns the closest plausible distribution, along with a boolean representing whether
    the original hypothesis was rejected. 
    """
    count_vector = [data.count(x) for x in value_list]
    s_prime, flipped = s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs)
    #the below print statements may need to be reworked
    if s_prime == 1:
        print("Proposed distribution not rejected at alpha = " + str(alpha) + ".")
        return (binary_hypothesis, False)
    elif flipped:
        q = [1-s_prime*binary_hypothesis[1], s_prime*binary_hypothesis[1]]
        print("Binomial tail exceeds 1 - "+ str(alpha) + ".")
    else:
        q = [s_prime*binary_hypothesis[0], 1-s_prime*binary_hypothesis[0]]
    q = [round(val, sigfigs-1) for val in q]
    print("Proposed distribution is rejected at alpha = " + str(alpha) + ". The closest plausible distribution is " + str(q) + ".")
    return (q, True)

#----------------------------------------------------------------------------------------------

#from mg_calculator_count_based import *

def kardis_test_main(counts, alpha, hypothesis):
    """
    Performs a multinomial kardis test upon the provided count vector 'counts'.
    'alpha' is the chosen alpha level, and 'hypothesis' is a list of 
    proposed probabilities for each value, 
    where counts[i] corresponds to hypothesis[i].
    Returns the kardis, a boolean corresponding to rejection, 
    as well as three values purely for use in the main function that calls this.
    """
    num_bins = len(counts)
    n = sum(counts)
    norm_scriptx = math.comb(n+num_bins-1, num_bins-1)
    mg = mg_calculator(counts, hypothesis)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    h = mpf(math.factorial(n)) / math.prod([mpf(math.factorial(x)) for x in counts])\
    * math.prod([mpf(hypothesis[x])**mpf(counts[x]) for x in range(num_bins)])
    kardis = r*h/nu
    if kardis < alpha:
        reject = True
    else:
        reject = False
    return (kardis, reject, r, nu, h)

def binary_kardis_test_main(counts, alpha, hypothesis):
    """
    Performs a binomial kardis test with tighter bounds than 'kardis_test_main'.
    All inputs are identical as in 'kardis_test_main', though 
    'counts' and 'hypothesis' may only have a length of 2.
    """
    n = sum(counts)
    k = counts[0]
    r = n+1
    p = mpf(hypothesis[0])**counts[0] * mpf(hypothesis[1])**counts[1]
    nu = mpf(math.comb(n, k))**(-1)
    kardis = r*p/nu
    if kardis < alpha:
        reject = True
    else:
        reject = False
    return (kardis, reject, r, nu, p)

def uniform_dist_kardis_test(observation, value_list, alpha):
    """
    Performs a multinomial kardis test using the uniform distribution.
    'observation' is the list of data points, 'value_list' is a list of the values
    a data point may take, and 'alpha' is the chosen alpha level. 
    Returns the kardis value and a boolean representing rejection.
    This function is largely deprecated.
    """
    obs_counts = [observation.count(x) for x in value_list]
    num_bins = len(value_list) #should == len(hyp) == len(obs_counts)
    #now |x| == num of counts 
    norm_scriptx = math.comb(len(observation)+num_bins-1, num_bins-1)
    u_hyp = [1/num_bins]*num_bins
    u = mpf(math.factorial(len(observation))) / math.prod([mpf(math.factorial(x)) for x in obs_counts])\
    * math.prod([mpf(u_hyp[x])**mpf(obs_counts[x]) for x in range(num_bins)])
    mg = mg_calculator(obs_counts, u_hyp)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    kardis = r*u/nu
    if kardis < alpha:
        reject = True
        print("Uniform distribution rejected at alpha = "+ str(alpha) + ". p(x) = " + str(u) + ". Kardis = " + str(kardis))
    else:
        reject = False
        print("Uniform distribution not rejected at alpha = "+ str(alpha)  + ". p(x) = " + str(u) +". Kardis = " + str(kardis))
    return (kardis, reject)

def univariate_kardis_test(observation, value_list, alpha, hypothesis=[]):
    """
    Performs a multinomial kardis test using a proposed distribution 'hypothesis'.
    If hypothesis is not entered, it will default to the uniform distribution. 
    'observation' is the list of data points, 'value_list' is a list of the values
    a data point may take, and 'alpha' is the chosen alpha level. 
    Returns the kardis value and a boolean representing rejection.
    This function is largely deprecated.
    """
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
    obs_counts = [observation.count(x) for x in value_list]
    num_bins = len(value_list) #should == len(hyp) == len(obs_counts)
    #now |x| == num of counts 
    norm_scriptx = math.comb(len(observation)+num_bins-1, num_bins-1)
    #don't forget the whole-number limitations on this calculator
    mg = mg_calculator(obs_counts, hyp)
    nu = norm_scriptx/mg
    r = norm_scriptx*(1+math.log(norm_scriptx))
    p = mpf(math.factorial(len(observation))) / math.prod([mpf(math.factorial(x)) for x in obs_counts])\
    * math.prod([mpf(hyp[x])**mpf(obs_counts[x]) for x in range(num_bins)])
    #print(p)
    kardis = r*p/nu
    if kardis < alpha:
        reject = True
        print("Proposed distribution rejected at alpha = "+ str(alpha)  + ". p(x) = " + str(p) +". Kardis = " + str(kardis))
    else:
        reject = False
        print("Proposed distribution not rejected at alpha = "+ str(alpha)  + ". p(x) = " + str(p) +". Kardis = " + str(kardis))
    return (kardis, reject)

#----------------------------------------------------------------------------------------------

def mg_calculator(observed_freq, hypothesis):
    """
    Calculates the number of count-based-vectors (events)
    that are more surprising than the observed event for
    any hypothesis. 'observed_freq' is the count vector 
    for the observation, and 'hypothesis' is the user-provided
    list of probabilities for each value. 
    
    Notes: The mean frequency of each value must be a whole number.
    That is, the each element of hypothesis multiplied by the observation
    length must be a whole number.
    """
    mean_freq = [int(sum(observed_freq)*i) for i in hypothesis]
    min_distance = sum(list(map(lambda x,y: abs(x-y), observed_freq, mean_freq)))
    max_distance_freq = [0 if i != mean_freq.index(min(mean_freq)) else sum(observed_freq) for i in range(len(observed_freq))]
    max_distance = sum(list(map(lambda x,y: abs(x-y), max_distance_freq, mean_freq)))
    num_bins = len(mean_freq)
    bin_list = bin_information(mean_freq)
    powerset_dict = {}
    for i in bin_list:
        powerset_dict[i[0]] = powerset_with_sums(i[0])
    mg = 0         
    for i in range(min_distance, max_distance+1, 2):
        if i == 0:
            mg += 1
        else:
            half_distance = i // 2
            valid_bins = list(filter(lambda x: x[2] >= half_distance, bin_list))
            for i in valid_bins:
                num_neg_bins = i[1]
                neg_placement_choices = num_sized_integer_compositions_multiple_limits(num_neg_bins, half_distance, powerset_dict[i[0]])        
                num_pos_bins = num_bins - num_neg_bins
                pos_placement_choices = num_weak_compositions(num_pos_bins, half_distance)
                mg += neg_placement_choices * pos_placement_choices
    return mg
                
def mg_calculator_uniform_hyp(observed_bin, mean):
    """
    Calculates the number of count-based-vectors (events)
    that are more surprising than the observed event for
    a uniform hypothesis
    
    Notes: The mean must be a whole number
    """
    min_distance = sum(list(map(lambda x: abs(x-mean), observed_bin)))
    max_distance = sum(observed_bin) - mean + (len(observed_bin)-1)*mean
    num_bins = len(observed_bin)
    mg = 0
    for i in range(min_distance, max_distance + 1, 2):
        if i == 0:
            mg += 1
        else:
            half_distance = i // 2
            min_neg_bins = math.ceil(half_distance/mean)
            max_neg_bins = num_bins - 1
            for num_neg_bins in range(min_neg_bins, max_neg_bins + 1):
                neg_bin_choices = math.comb(num_bins, num_neg_bins)
                neg_placement_choices = num_sized_integer_compositions_uniform_limit(num_neg_bins, half_distance, mean)
                num_pos_bins = num_bins - num_neg_bins
                pos_placement_choices = num_weak_compositions(num_pos_bins, half_distance)
                mg += neg_bin_choices * neg_placement_choices * pos_placement_choices
    return mg

def num_weak_compositions(length, total):
    """
    Calculates how many ways there are to distribute n
    balls into k bins (allowing for empty bins)
    
    See: https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics)
    """
    k, N = length, total
    return math.comb(N+k-1, k-1)
    
def num_sized_integer_compositions_uniform_limit(length, total, limit):
    """
    Calculates how many ways there are to distribute N
    balls into n bins (not allowing for empty bins)
    where each bin can hold at maximum r balls
    
    See: https://math.stackexchange.com/questions/553960/extended-stars-and-bars-problemwhere-the-upper-limit-of-the-variable-is-bounded
    """
    n, N, r = length, total, limit
    end = min(length, total//(limit+1))
    return sum([((-1)**k)*(math.comb(n, k))*(math.comb(N-k*(r)-1, n-1)) for k in range(end+1)])

def num_sized_integer_compositions_multiple_limits(length, total, powerset_list):
    """
    Calculates how many ways there are to distribute N
    balls into n bins (not allowing for empty bins)
    where each bin has a unique maximum number of balls
    it can hold

    See: https://math.stackexchange.com/questions/553960/extended-stars-and-bars-problemwhere-the-upper-limit-of-the-variable-is-bounded
    """
    n, N = length, total
    composition_count = 0
    for i in powerset_list:
        m = N-1-i[2]
        if m >= 0:
            composition_count += (-1)**i[1] * math.comb(m, n-1)
    return composition_count

def bin_information(bin_list):
    all_possible_bin_combos = []
    for i in range(len(bin_list)):
        all_possible_bin_combos += list(combinations(bin_list, i))
    bin_list = [(i, len(i), sum(i)) for i in all_possible_bin_combos]
    return bin_list

def powerset_with_sums(orig_set):
    powerset = []
    for i in range(len(orig_set)+1):
        powerset += list(combinations(orig_set, i))
    powerset_with_sums = [[i, len(i), sum(i)] for i in powerset]
    return powerset_with_sums

#----------------------------------------------------------------------------------------------

# Used to suppress useless Scipy warnings
# Disable this if you're rigorously testing
warnings.filterwarnings("ignore")

def q_finder_main(count_vector, hypothesis, p_lowerbound):
    """
    Returns the closest plausible distribution to the 'hypothesis' distribution. 
    'count_vector' contains the quantity of each value in the data. 
    'p_lowerbound' is the lower bound on plausible probabilities - any
    distribution which produces a lower sequence probability
    may be instantly rejected.
    """
    p = np.array(hypothesis)
    log_n_fac = math.log(math.factorial(sum(count_vector)))
    log_prod_x_fac = math.log(math.prod([math.factorial(x) for x in count_vector]))
    constraint_constant = math.log(p_lowerbound) - log_n_fac + log_prod_x_fac
    # Defines loss function to be the KL divergence of Q from P
    def objective(q):
        return sum(rel_entr(q, p))
    # The sum of all probabilities in Q must be 1
    linear_constraint = LinearConstraint(len(count_vector)*[1], [1], [1], keep_feasible=False)
    # log(Q(x)) - constraint_constant >= 0
    def nonlin_cons(q):
        log_prod_q = sum([count_vector[i]*math.log(q[i]) for i in range(len(count_vector))])
        constraint = log_prod_q - constraint_constant
        return constraint
    nonlinear_constraint = NonlinearConstraint(nonlin_cons, lb=0, ub=np.inf, keep_feasible=True)
    # Each probability in Q must be between 0 and 1 inclusive
    bounds = Bounds(len(count_vector)*[0], len(count_vector)*[1], keep_feasible=True)
    # Bogus initial guess for the algorithm to start with
    x0 = np.array(len(count_vector)*[0.9])
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint, nonlinear_constraint], 
                   options={'gtol': 1e-100, 'xtol': 1e-100, 'maxiter': 1000000, 'verbose': 0}, bounds=bounds)
    return res.x

def binary_q_finder_main(binary_count_vector, hypothesis, alpha):
    """
    Acts similarly to 'q_finder_main', but is strictly for binary situations.
    'hypothesis' is the originally proposed distribution, 
    and 'alpha' is the given alpha level. 
    """
    p = np.array(hypothesis)
    n = sum(binary_count_vector)
    k = binary_count_vector[0]
    constraint_constant = math.log(alpha) - math.log(n+1) - math.log(math.comb(n,k))
    def objective(q):
        return sum(rel_entr(q, p))
    linear_constraint = LinearConstraint([1,1], [1], [1], keep_feasible=False)
    def nonlin_cons(q):
        return k*math.log(q[0]) + (n-k)*math.log(q[1]) - constraint_constant
    nonlinear_constraint = NonlinearConstraint(nonlin_cons, lb=0, ub=np.inf, keep_feasible=True)
    bounds = Bounds([0,0], [1,1], keep_feasible=True)
    x0 = np.array([0.99,0.99])
    res = minimize(objective, x0, method='trust-constr', constraints=[linear_constraint, nonlinear_constraint], 
                   options={'gtol': 1e-100, 'xtol': 1e-100, 'maxiter': 1000000, 'verbose': 0}, bounds=bounds)
    return res.x

#----------------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------------

#from data_binarizer import *


def s_prime_finder_main(count_vector, value_list, selected_value_list, alpha, binary_hypothesis, sigfigs):
    """
    For the tightest-bound binary function. Finds the coefficient necessary 
    on the selected value probability to produce a plausible distribution. 
    'count_vector' is the original non-binary count representation of the data,
    'value_list' contains all the values a data point may take, 
    'selected_value_list' is the collection of values the user believes
    the data is biased towards. 'alpha' is the requested alpha level, 
    'binary_hypothesis' contains the probability of producing a
    selected value or a not-selected value,
    and 'sigfigs' is the requested number of significant figures.
    """
    binom_dict = {}
    flipped = False
    bcounts = data_binarizer_main(count_vector, value_list, selected_value_list)
    n = sum(count_vector)
    i = bcounts[0]
    p_b = binary_hypothesis[0]

    def stirl_approx(k):
            """
            Approximates nCk using Stirling's approximation.
            """
            return sqrt(n)*mpf(n)**mpf(n)/(sqrt(2*math.pi*k*(n-k))*mpf(k)**k*mpf(n-k)**(n-k))
    def binom(p):
        """
        Returns the binomial one-sided (greater) tail probability.
        """
        if p in binom_dict.keys():
            return binom_dict[p]
        else:
            prob_more_extreme = sum([stirl_approx(k)*(mpf(p)**k)*(mpf(1-p)**(n-k)) for k in range(i,n)])
            prob_more_extreme += mpf(p)**n
            binom_dict[p] = prob_more_extreme
            return prob_more_extreme

    if binom(p_b) >= alpha:
        if binom(p_b) > 1-alpha:
            flipped = True
            i = bcounts[1]
            p_b = binary_hypothesis[1]
        else:
            return 1, flipped
    cur_pow10 = 0
    cur_val = round(1/p_b, sigfigs-1)

    def buddy(s):
        """
        Recursively optimizes s' for a single power of 10.
        """
        lo = (s - 10**cur_pow10)
        lo_binom = binom(lo*p_b)
        if lo_binom < alpha:
            return s
        else:
            return buddy(lo)
            
    for x in range(sigfigs):
        cur_val = buddy(cur_val)
        cur_pow10 -= 1

    return (round(cur_val,sigfigs-1),flipped)

#----------------------------------------------------------------------------------------------

def data_slicer_multival(nbdata, nbhypothesis, nbvaluelist, selectedvaluelist):
    
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

#----------------------------------------------------------------------------------------------

def graph_times(num_bins, avg_run_times):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.font_manager as font_manager

    plt.figure(dpi=300)

    plt.plot(num_bins,avg_run_times,'-r*')
    plt.xlabel("Number of Bins",fontsize = 14,fontname="Sans-serif")
    plt.ylabel("Time (seconds)",fontsize = 14,fontname="Sans-serif")
    plt.xticks(num_bins)
    plt.title("Hypothesis Test Average Run Time per Number of Bins",fontsize = 14,fontname="Sans-serif")
    plt.ylim([0,max(avg_run_times)+.5])

    #legend_font = font_manager.FontProperties(family='sans-serif',style='normal', size=10)
    #plt.legend(prop=legend_font)
    plt.grid(True)
    plt.savefig("Time scaling experiment results graph.pdf",bbox_inches="tight")
    
#----------------------------------------------------------------------------------------------