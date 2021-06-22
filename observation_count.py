"""
Based off code from https://codereview.stackexchange.com/questions/202518/count-occurrences-of-a-specific-sequence-in-a-list-of-many-numbers

This function returns the count vector(list of counts) for the observed sequence

The observation and value_list should be in the format of a list of list
i.e for a 6-sided die:  
observation = 45*[[1]]+5*[[2]]+4*[[3]]+3*[[4]]+2*[[5]]+[[6]] 
value_list = [[1],[2],[3],[4],[5],[6]]

Output: [45, 5, 4, 3, 2, 1]

This formating rule also applies to multivariate observations and value_lists.
"""


def observation_count(observation, value_list):
    flatten_obs = [i for x in observation for i in x]
    count_vector = []
    for value in value_list:
        count = 0
        num_classes = len(value)
        upper_bound = len(flatten_obs)-num_classes+1
        for i in range(upper_bound):
             if flatten_obs[i:i+num_classes] == value:
                count += 1
        count_vector.append(count)
    return count_vector
