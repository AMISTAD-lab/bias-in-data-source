import math
from itertools import *

def mg_calculator(observed_freq, hypothesis):
    """
    Calculates the number of count-based-vectors (events)
    that are more surprising than the observed event for
    any hypothesis. 'observed_freq' is the count vector 
    for the observation, and 'hypothesis' is the user-provided
    list of probabilities for each value.
    """
    mean_freq = [int(round(sum(observed_freq)*i, 0)) for i in hypothesis]
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
                neg_placement_choices = num_sized_integer_compositions(num_neg_bins, half_distance, powerset_dict[i[0]])        
                num_pos_bins = num_bins - num_neg_bins
                pos_placement_choices = num_weak_compositions(num_pos_bins, half_distance)
                mg += neg_placement_choices * pos_placement_choices
    return mg

def num_weak_compositions(length, total):
    """
    Calculates how many ways there are to distribute n
    balls into k bins (allowing for empty bins)
    
    See: https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics)

    This is a Lemma in the paper.
    """
    k, N = length, total
    return math.comb(N+k-1, k-1)

def num_sized_integer_compositions(length, total, powerset_list):
    """
    Calculates how many ways there are to distribute N
    balls into n bins (not allowing for empty bins)
    where each bin has a unique maximum number of balls
    it can hold

    This is a Theorem in the paper. 
    """
    n, N = length, total
    composition_count = 0
    for i in powerset_list:
        m = N-1-i[2]
        if m >= 0:
            composition_count += (-1)**i[1] * math.comb(m, n-1)
    return composition_count

def bin_information(bin_list):
    """
    For a given list of possible bins where each bin is simply represented
    as a non-negative integer corresponding to its size, this function returns
    a list of bin combination tuples. Each tuple contains a combination of bins,
    the number of bins in that combination, and the total capacity of the bins.
    """
    all_possible_bin_combos = []
    for i in range(len(bin_list)):
        all_possible_bin_combos += list(combinations(bin_list, i))
    bin_list = [(i, len(i), sum(i)) for i in all_possible_bin_combos]
    return bin_list

def powerset_with_sums(orig_set):
    """
    For a given set (combination) of bins (defined identially in `bin_information`),
    this function returns a list of lists where each sublist contains an element of
    the original set's powerset, that elements length, and total capacity.

    This function is used to precompute the powerset of all bin combinations to reduce
    the computational runtime of `mg_calculator`.
    """
    powerset = []
    for i in range(len(orig_set)+1):
        powerset += list(combinations(orig_set, i))
    powerset_with_sums = [[i, len(i), sum(i)] for i in powerset]
    return powerset_with_sums