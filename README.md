# bias-in-data-source
Repo for the ```Identifying Bias in Data using Two-Distribution Hypothesis Tests``` article.

What's here:

## hypothesis_test.py
This file contains several functions to help you run a two-distribution specified complexity (s.c.) hypothesis test for a discrete univariate distribution. This is the user-facing file, and we recommend using it instead of the other files in this directory. These functions share many of the same inputs, which will only be described once. 

#### hypothesis_test(data, value_list, alpha = 0.05, hypothesis = [])
Performs a multinomial SC test.
  - data: a list representing a sequence of events observed
  - value_list: a list of all possible values the discrete variable in question can take on
  - alpha: the alpha level for your hypothesis test. The default is 0.05.
  - hypothesis: a list representing a probability distribution for each of the values in value_list. The default is uniform.
  
#### Use example
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

## Other Files
The other files in the main folder can be split into two categories:
#### Support for hypothesis_test.py
  - kardis_test.py
  - mg_calculator.py
  - q_finder.py
#### Experiments and data visualization
These ipynb files can be run to recreate the tests we performed. 
  - compas_experiments_recid_score.ipynb
  - compas_experiments_charge_type.ipynb
  - uci_adult_experiments.ipynb
  - run_time_experiments_decile_score.ipynb
  - run_time_experiments_score_text.ipynb
  - graph_distributions.py
