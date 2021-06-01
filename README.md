# Bias-in-Data
Repo for the Bias in Data summer project

What's here so far:

## univariate_sc_test.py
This file contains several functions to help you run a two-distribution specified complexity (s.c.) hypothesis test for a discrete univariate distribution. The main function that allows you to do this is univariate_sc_test(observation, value_list, hypothesis, alpha).
```
observation - a list representing a sequence of events observed (see example below)
value_list - a list of all possible values the discrete variable in question can take on
hypothesis - a list representing a probability distribution for each of the values in value list
alpha - the alpha level for your hypothesis test
```
### How to use
If you flipped a coin (which you hypothesized was fair) 20 times and got a sequence of 20 heads in a row, you could conduct a s.c. test to determine whether your hypothesis could be a plausible explanation for the event you observed. If you were to pick an alpha of 0.05, the parameters representing this scenario would be 
```
observation = 20*['H']
value_list = ['H', 'T']
hypothesis = [0.5, 0.5]
alpha = 0.05
```
and running the function would look like
```
univariate_sc_test(20*['H'], ['H', 'T'], [0.5, 0.5], 0.05]
```
The function returns a tuple contain p(x), s\*u(x), and boolean representing whether the hypothesis was rejected or not. 
If the hypothesis (proposed distribution) is rejected by the s.c. test, the function will also print 
