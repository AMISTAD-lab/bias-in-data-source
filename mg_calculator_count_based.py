import math
from itertools import *

def mg_calculator(observed_freq, hypothesis):
    """
    Calculates the number of count-based-vectors (events)
    that are more surprising than the observed event for
    any hypothesis
    
    Notes: The mean frequency of each value must be a whole number.
    That is, the each element of hypothesis multiplied by the observation
    length must be a whole number.
    """
    mean_freq = [int(sum(observed_freq)*i) for i in hypothesis]
    min_distance = sum(list(map(lambda x,y: abs(x-y), observed_freq, mean_freq)))
    max_distance_freq = [0 if i != mean_freq.index(min(mean_freq)) else sum(observed_freq) for i in range(len(observed_freq))]
    max_distance = sum(list(map(lambda x,y: abs(x-y), max_distance_freq, mean_freq)))
    num_bins = len(mean_freq)
    bin_dict, bin_list = bin_information(mean_freq)
    bin_list.sort(key=lambda x: x[2])
    powerset_dict = {}
    for i in bin_list:
        powerset_dict[i[0]] = powerset_with_sums(i[0])
    mg = 0         
    for i in range(min_distance, max_distance+1, 2):
        half_distance = i // 2
        valid_bins = list(filter(lambda x: x[2] >= half_distance, bin_list))
        for i in valid_bins:
            limit_list = i[0]
            num_neg_bins = i[1]
            neg_placement_choices = num_sized_integer_compositions_multiple_limits(num_neg_bins, half_distance, limit_list, powerset_dict[i[0]])        
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

def num_sized_integer_compositions_multiple_limits(length, total, limit_list, powerset_list):
    """
    Calculates how many ways there are to distribute N
    balls into n bins (not allowing for empty bins)
    where each bin has a unique maximum number of balls
    it can hold

    See: https://math.stackexchange.com/questions/553960/extended-stars-and-bars-problemwhere-the-upper-limit-of-the-variable-is-bounded
    """
    n, N, r = length, total, limit_list
    composition_count = 0
    for i in powerset_list:
        m = N-1-i[2]
        if m > 0:
            composition_count += (-1)**i[1] * math.comb(m, n-1)
    return composition_count

def bin_information(bin_list):
    all_possible_bin_combos = []
    for i in range(len(bin_list)):
        all_possible_bin_combos += list(combinations(bin_list, i))
    bin_dict = {}
    for i in all_possible_bin_combos:
        bin_dict[i] = sum(i)
    bin_list = [(i, len(i), sum(i)) for i in all_possible_bin_combos]
    return bin_dict, bin_list

def powerset_with_sums(orig_set):
    powerset = []
    for i in range(len(orig_set)+1):
        powerset += list(combinations(orig_set, i))
    powerset_with_sums = [(i, len(i), sum(i)) for i in powerset]
    return powerset_with_sums