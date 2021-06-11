import math
from itertools import *

"""
The sized_partitions function is from
https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning
"""

def fast_mg_calculator(observation, value_list, hypothesis):
    """
    Calculates Mg(x) for a given observation and list of 
    possible values for a given discrete random variable
    utilizing sum of absolute distances as a difference 
    measure between two distributions
    """
    obs_bin = [observation.count(i) for i in value_list]
    mean_bin = list(map(lambda x: x*len(observation), hypothesis))
    min_distance = math.fsum(list(map(lambda x,y: abs(x-y), obs_bin, mean_bin)))
    mg = 0
    scriptx = scriptx_generator(len(value_list), len(observation))
    event_list = []
    for event in scriptx:
        test_distance = math.fsum(list(map(lambda x,y: abs(x-y), event, mean_bin)))
        if test_distance >= min_distance:
            mg += math.factorial(len(observation)) // math.prod(list(map(lambda x: math.factorial(x), event)))
            event_list.append(event)
    return mg, event_list

def scriptx_generator(num_vals, observation_length):
    return permute_events(subtract_1_from_all(list(sized_partitions(observation_length+num_vals, num_vals))))

def sized_partitions(n, k, m = None):
  """Partition n into k parts with a max part of m.
  Yield non-increasing lists.  m not needed to create generator.
  """
  if k == 1:
    yield [n]
    return
  for f in range(n-k+1 if (m is None or m > n-k+1) else m, (n-1)//k, -1): 
    for p in sized_partitions(n-f, k-1, f): yield [f] + p
    
def subtract_1_from_all(l):
    for i in range(len(l)):
        l[i] = list(map(lambda x: x-1, l[i]))
    return l

def permute_events(l):
    all = []
    for i in range(len(l)):
        all += list(set(list(permutations(l[i]))))
    return all