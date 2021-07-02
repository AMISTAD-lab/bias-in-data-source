# Bias-in-Data
Repo for the Bias in Data summer project

What's here so far:

## hypothesis_test.py
This file contains several functions to help you run a two-distribution specified complexity (s.c.) hypothesis test for a discrete univariate distribution. This is the user-facing file, and we recommend using it instead of the other files in this directory.

#### hypothesis_test(data, value_list, alpha = 0.05, hypothesis = [])
Performs a multinomial SC test.
##### How to use:

#### binary_hypothesis_test(data, value_list, selected_value_list, alpha=0.05, binary_hypothesis=[0.5,0.5])
Performs a binomial SC test with tighter bounds than hypothesis_test.
##### How to use:

#### exact_binomial_test(data, value_list, selected_value_list, alpha = 0.05, binary_hypothesis = [0.5,0.5], sigfigs = 4)
Performs an exact binomial test. Has the tightest bounds of all the functions listed, but is significantly slower.
##### How to use:
