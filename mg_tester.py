import math
from scipy.special import rel_entr

def mg_tester(observation, value_list, hypothesis):
    """
    Temporary function for testing Mg(x)
    """
    obs_dist = [observation.count(i)/len(observation) for i in value_list]
    min_kl = math.fsum(rel_entr(obs_dist, hypothesis))
    mg = 0
    hyp_bin = list(map(lambda x: x*len(observation), hypothesis))
    obs_bin = list(map(lambda x: x*len(observation), obs_dist))
    orig_distance = list(map(lambda x,y: x-y, obs_bin, hyp_bin))
    orig_sum_of_squares = math.fsum([i**2 for i in orig_distance])
    max_orig_dist = max(orig_distance)
    print(str(orig_distance) + " : " + str(orig_sum_of_squares))
    scriptx = scriptx_generator_helper(scriptx_generator(len(value_list), len(observation)), len(value_list))
    for event in scriptx:
        test_dist = [event[i]/len(observation) for i in range(len(value_list))]
        test_kl = math.fsum(rel_entr(test_dist, hypothesis))
        distance = list(map(lambda x,y: x-y, event, hyp_bin))
        abs_dist = [abs(i) for i in distance]
        test_sum_of_squares = math.fsum([i**2 for i in distance])
        if test_kl >= min_kl:
            """
            if test_sum_of_squares < orig_sum_of_squares:
                print(str(event) + ' : ' + str(distance) + ' : ' + str(test_sum_of_squares))
            mg += math.factorial(len(observation)) // math.prod(list(map(lambda x: math.factorial(x), event)))
            """
        else:
            if 0 in event:
                print(str(event) + ' : ' + str(abs_dist))
    return mg

def scriptx_generator(num_vals, observation_length, current_vals=[]):
    """
    Generates a representation of scriptx
    for a discrete random variable
    Notes:
        This representation of scriptx is based
        "bins" of values, where scriptx is a list
        of lists. Each list in scriptx is a
        possible frequency distribution. However,
        this function itself does not do the grouping
        of these distributions into lists. Rather, this
        function outputs just a single list of numbers
        that are grouped into distributions by the helper
        function below.
    Example: Rolling a fair die 6 times
        Some of the possible frequency distributions
        in scriptx are [1,1,1,1,1,1], [0,6,0,0,0,0],
        and [1,0,3,2,0,0]. These would show up in scriptx
        as the single list 
        [1,1,1,1,1,1,0,6,0,0,0,0,1,0,3,2,0,0]
        and would be grouped into nice lists by the helper
        function below.
    """
    if observation_length == 0:
        return current_vals + [0]*num_vals
    elif num_vals == 1:
        return current_vals + [observation_length]
    else:
        scriptx = []
        for i in range(observation_length+1):
            scriptx += scriptx_generator(num_vals-1, observation_length-i, current_vals + [i])
        return scriptx

def scriptx_generator_helper(ungrouped_scriptx_list, num_vals):
    """
    Helper function for scriptx_generator
    Notes:
        scriptx generator only generates a single
        list of numbers which are in the correct order,
        but still need to be grouped into frequency
        distributions in order to be interpretable.
        This function does this grouping.
    """
    grouped_perm_list = []
    for i in range(0, len(ungrouped_scriptx_list), num_vals):
        grouped_perm_list.append(ungrouped_scriptx_list[i:i+num_vals])
    return grouped_perm_list