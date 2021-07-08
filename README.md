# Bias-in-Data
Repo for the Bias in Data summer project

What's here so far:

## hypothesis_test.py
This file contains several functions to help you run a two-distribution specified complexity (s.c.) hypothesis test for a discrete univariate distribution. This is the user-facing file, and we recommend using it instead of the other files in this directory. These functions share many of the same inputs, which will only be described once. 

#### hypothesis_test(data, value_list, alpha = 0.05, hypothesis = [])
Performs a multinomial SC test.
  - data: a list representing a sequence of events observed
  - value_list: a list of all possible values the discrete variable in question can take on
  - alpha: the alpha level for your hypothesis test. The default is 0.05.
  - hypothesis: a list representing a probability distribution for each of the values in value_list. The default is uniform.

#### binary_hypothesis_test(data, value_list, selected_value_list, alpha=0.05, binary_hypothesis=[0.5,0.5])
Performs a binomial SC test with tighter bounds than hypothesis_test. 
  - selected_value_list: a list of values that is used to binarize the data, which will be partitioned into selected values and not-selected values. 
  - binary_hypothesis: a list representing a binary probability distribution for selected values and not-selected values. The default is uniform.

#### exact_binomial_test(data, value_list, selected_value_list, alpha = 0.05, binary_hypothesis = [0.5,0.5], sigfigs = 4)
Performs an exact binomial test. Has the tightest bounds of all the functions listed, but is significantly slower. Like the above binary test, the provided data should still be non-binary. 
  - sigfigs: the number of significant figures the returned distribution will have. The default is 4. 
  
#### Use examples
Suppose you have a six-sided die which you believed to be fair. You rolled it sixty times and observed 45 ones, 5 twos, 4 threes, 3 fours, 2 fives, and 1 six. If you were to pick an alpha of 0.05, the parameters representing this scenario would be 
  - data: 45*[1] + 5*[2] + 4*[3] + 3*[4] + 2*[5] + [6]
  - value_list: [1, 2, 3, 4, 5, 6]
  - alpha: 0.05
  - hypothesis: 6*[1/6] 
Running the hypothesis_test function would look like 
```
hypothesis_test(45*[1]+5*[2]+4*[3]+3*[4]+2*[5]+[6], [1,2,3,4,5,6], 0.05, 6*[1/6])
```
or, more simply, 
```
hypothesis_test(45*[1]+5*[2]+4*[3]+3*[4]+2*[5]+[6], [1,2,3,4,5,6])
```
as the chosen alpha and hypothesis values are the default. This function will return the kardis value and a boolean representing whether the hypothesis was rejected or not.
If the hypothesis was rejected, then the function will additionally return the probability coefficient needed to produce a plausible explanation, the lower plausibility bound that results from that coefficient, and the closest plausible distribution. For the die example above, the SC test would determine that the uniform hypothesis was not a reasonable explanation, and the printed distribution would be (here rounded to three significant figures) 
```
[0.466, 0.133, 0.121, 0.108, 0.094, 0.078]
```
which indicates the die would need to have at least a 46.6% chance of rolling a 1. 

If you instead wanted to use our other tests, assuming that the die was biased towards 1, the parameters would be 
  - data: 45*[1] + 5*[2] + 4*[3] + 3*[4] + 2*[5] + [6]
  - value_list: [1, 2, 3, 4, 5, 6]
  - selected_value_list: [1]
  - alpha: 0.05
  - binary_hypothesis[1/6, 5/6]

```
binary_hypothesis_test(45*[1]+5*[2]+4*[3]+3*[4]+2*[5]+[6], [1,2,3,4,5,6], [1], 0.05, [1/6,5/6])
```
Returns a closest plausible distribution of ```[0.553, 0.447]```, which means a 55.3% chance of rolling a 1 is required. 
```
exact_binomial_test(45*[1]+5*[2]+4*[3]+3*[4]+2*[5]+[6], [1,2,3,4,5,6], [1], 0.05, [1/6,5/6])
```
Returns a closest plausible distribution ```[0.641, 0.359]```, which means a 64.1% chance of rolling a 1 is required.
This illustrates the slightly different bounds on each test. 

## Other Files
The other files in the main folder can be split into two categories:
#### Support for hypothesis_test.py
  - counts_kardis_test.py
  - data_binarizer.py
  - mg_calculator_count_based.py
  - q_finder_count_based.py
  - s_prime.py
#### Data processing and visualization
  - csv_to_lists.py
  - data_bin_slicer_multival.py
  - excel_to_list.py
  - graph_distributions.py

