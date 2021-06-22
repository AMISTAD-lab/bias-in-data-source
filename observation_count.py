"""
Based off code from https://codereview.stackexchange.com/questions/202518/count-occurrences-of-a-specific-sequence-in-a-list-of-many-numbers

This function returns the count vector for the observed sequence

value_list should be in the format of a list of list i.e for a 6-sided
die the value_list will be  [[1],[2],[3],[4],[5],[6]]
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
