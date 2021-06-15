import math
from itertools import *

def mg_calculator_event_based(observed_freq, hypothesis):
    """
    Calculates the number of count-based-vectors (events)
    that are more surprising than the observed event for
    any hypothesis
    
    Notes: The mean frequency of each value must be a whole number.
    That is, the each element of hypothesis multiplied by the observation
    length must be a whole number.
    """
    mean_freq = [int(sum(observed_freq)*i) for i in hypothesis]
    indexed_mean_freq = [(i, mean_freq[i]) for i in range(len(mean_freq))]
    min_distance = sum(list(map(lambda x,y: abs(x-y), observed_freq, mean_freq)))
    max_distance_freq = [0 if i != mean_freq.index(min(mean_freq)) else sum(observed_freq) for i in range(len(observed_freq))]
    max_distance = sum(list(map(lambda x,y: abs(x-y), max_distance_freq, mean_freq)))
    num_bins = len(observed_freq)
    event_count = 0
    for i in range(min_distance, max_distance+1, 2):
        half_distance = i // 2
        min_neg_bins = min_bins_required(mean_freq, half_distance)
        max_neg_bins = num_bins - 1
        for num_neg_bins in range(min_neg_bins, max_neg_bins+1):
            neg_bin_choices = feasible_bin_choices(indexed_mean_freq, num_neg_bins, half_distance)
            for neg_bins in neg_bin_choices:
                limit_list = [i[1] for i in neg_bins]
                neg_placement_choices = num_partitions_multiple_limits(len(neg_bins), half_distance, limit_list)
                num_pos_bins = len(indexed_mean_freq) - len(neg_bins)
                pos_placement_choices = num_partitions_unconstrained(num_pos_bins, half_distance)
                event_count += neg_placement_choices * pos_placement_choices
    return event_count
                
def mg_calculator_event_based_uniform_hyp(observed_bin, mean):
    """
    Calculates the number of count-based-vectors (events)
    that are more surprising than the observed event for
    a uniform hypothesis
    
    Notes: The mean must be a whole number
    """
    min_distance = sum(list(map(lambda x: abs(x-mean), observed_bin)))
    max_distance = sum(observed_bin) - mean + (len(observed_bin)-1)*mean
    num_bins = len(observed_bin)
    event_count = 0
    for i in range(min_distance, max_distance + 1, 2):
        half_distance = i // 2
        min_neg_bins = math.ceil(half_distance/mean)
        max_neg_bins = num_bins - 1
        for num_neg_bins in range(min_neg_bins, max_neg_bins + 1):
            neg_bin_choices = math.comb(num_bins, num_neg_bins)
            neg_placement_choices = num_partitions_uniform_limit(num_neg_bins, half_distance, mean)
            num_pos_bins = num_bins - num_neg_bins
            pos_placement_choices = num_partitions_unconstrained(num_pos_bins, half_distance)
            event_count += neg_bin_choices * neg_placement_choices * pos_placement_choices
    return event_count

def num_partitions_unconstrained(partition_size, total):
    """
    Calculates how many ways there are to distribute n
    balls into k bins (allowing for empty bins)
    
    See: https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics)
    """
    k, n = partition_size, total
    return math.comb(n+k-1, k-1)
    
def num_partitions_uniform_limit(partition_size, total, limit):
    """
    Calculates how many ways there are to distribute N
    balls into n bins (not allowing for empty bins)
    where each bin can hold at maximum r balls
    
    See: https://math.stackexchange.com/questions/553960/extended-stars-and-bars-problemwhere-the-upper-limit-of-the-variable-is-bounded
    """
    n, N, r = partition_size, total, limit
    end = min(partition_size, total//(limit+1))
    return sum([((-1)**q)*(math.comb(n, q))*(math.comb(N-q*(r)-1, n-1)) for q in range(end+1)])

def num_partitions_multiple_limits(partition_size, total, limit_list):
    """
    Calculates how many ways there are to distribute N
    balls into n bins (not allowing for empty bins)
    where each bin has a unique maximum number of balls
    it can hold

    See: https://math.stackexchange.com/questions/553960/extended-stars-and-bars-problemwhere-the-upper-limit-of-the-variable-is-bounded
    """
    n, N, r = partition_size, total, limit_list
    powerset_list = powerset(list(range(1, n+1)))
    partition_count = 0
    for i in powerset_list:
        m = N-1-sum([r[j-1] for j in i])
        k = n-1
        if m >= 0:
            partition_count += (-1)**len(i) * math.comb(m, k)
    return partition_count

def powerset(original_set):
    """
    Generates the powerset of original_set,
    which is to be provided as a list
    """
    powerset_list = []
    for i in range(len(original_set)+1):
        powerset_list += list(combinations(original_set, i))
    return powerset_list

def min_bins_required(bin_list, lower_limit):
    """
    Returns the minimum number of bins from a list 
    that are required to hold lower_limit balls
    """
    if max(bin_list) >= lower_limit:
        return 1
    else:
        max_val = max(bin_list)
        temp_bin_list = bin_list[:]
        temp_bin_list.remove(max_val)
        return 1 + min_bins_required(temp_bin_list, lower_limit-max_val)
    
def feasible_bin_choices(indexed_bin_list, num_bins, lower_limit):
    """
    Returns a list of all feasible combinations of num_bin
    bins that can hold lower_limit balls
    """
    all_bin_choices = list(combinations(indexed_bin_list, num_bins))
    return list(filter(lambda x: sum([i[1] for i in x]) >= lower_limit, all_bin_choices))