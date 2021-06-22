# Bias-in-Data
Repo for the Bias in Data summer project

What's here so far:

## counts_sc_test.py
This file contains several functions to help you run a two-distribution specified complexity (s.c.) hypothesis test for a discrete univariate distribution. The main function that allows you to do this is univariate_sc_test(observation, value_list, hypothesis, alpha).
- observation: a list representing a sequence of events observed (see example below)
- value_list: a list of all possible values the discrete variable in question can take on
- hypothesis: a list representing a probability distribution for each of the values in value list
- alpha: the alpha level for your hypothesis test

### How to use
If you had an electronic coin flipper (which you hypothesized was fair) that flipped a coin 20 times and landed a sequence of 20 heads in a row, you could conduct a s.c. test to determine whether your hypothesis could be a plausible explanation for the event you observed. If you were to pick an alpha of 0.05, the parameters representing this scenario would be 
- observation = 20*['H']
- value_list = ['H', 'T']
- hypothesis = [0.5, 0.5]
- alpha = 0.05
and running the function would look like
```
univariate_sc_test(20*['H'], ['H', 'T'], [0.5, 0.5], 0.05]
```
The function returns a tuple containing p(x), s\*u(x), and boolean representing whether the hypothesis was rejected or not.

If the hypothesis (proposed distribution) is rejected by the s.c. test, the function will additionally print another distribution, representing the closest distribution to the hypothesis which also does not get rejected by the s.c. test. In essence, this distribution puts a lower bound on the distribution such that it will pass the s.c. test. For the coin example above, the s.c. test would determine that the fair flipper hypothesis could not be a reasonble explanation, and the printed distribution would be
```
[0.7266, 0.2734]
```
which indicates that the coin flipper would have to have at least a 72.66% chance of landing heads in order to be considered as a plausible explanation for the event of 20 heads in a row.
