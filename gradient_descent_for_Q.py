import itertools
from sklearn.linear_model import LogisticRegression
import numpy as np
from numpy.random import default_rng

"""Gender is represented by the int values: 1 = M, 2 = F, 3 = Nb, 4 = Other. 
   Race is represented by the int values: 4 = White, 5 = Black, 6 = Asian,
   7 = Latino, 8 = Other"""
"""
Below is what I used to create a test dataset of 20 people with a specific gender and race.

rng = default_rng()
variablelist = [1, 2, 3, 4]
probabilitylist = [0.7, 0.2, 0.05, 0.05] #0.7 corresponds to m, 0.2 to f, etc
size = 20
datasetg = rng.choice(variablelist, size = size, p = probabilitylist) 

variablelist = [4, 5, 6, 7, 8]
probabilitylist = [0.5, 0.1, 0.3, 0.08, 0.02] #0.7 corresponds to m, 0.2 to f, etc
size = 20
datasetr = rng.choice(variablelist, size = size, p = probabilitylist) 
data = np.array(list(zip(datasetg,datasetr)))

#This below randomizes whether each person was hired or not

#variablelist = [0, 1]
#probabilitylist = [0.7, 0.3]
#size = 20
#outcomes = rng.choice(variablelist, size = size, p = probabilitylist) 
outcomes = []

#Below is used to only have white males be hired

for x in data:
    if x[0] == 1 and x[1] == 4:
        outcomes.append(1)
    else:
        outcomes.append(0)
"""

"""Function requires the the complete data set, a corresponding set of the class outcomes 
   for the raw data (i.e For hiring data, give a set saying if each person was hired or not), 
   and lastly a dictionary of the variables taken into account for each data point and 
   the numerical values is needed. """

def gradient_descent_for_Q(data,outcomes,variables):
    X = data
    y = outcomes
    clf = LogisticRegression(random_state=0)
    clf.fit(X,y)
    
    var = variables
    """Turns each key/value pair in the dictionary into separate arrays"""
    for key in var.keys():
        x = key+' = '+str(var[key])
        exec(x,globals())
        
    """Creates all the possible combinations of the variables. The variable names 
       placed in prodcut should be the names of all the keys in the variables dictionary. 
       Currently the default is set to gender and race since the example used had those as 
       the key names in variables."""  
    p = list(itertools.product(gender,race))
    combinations = np.array(p)
    
    """Gives an array of the probability of each combinations being placed in each class(i.e for our
    examples its whether they are hired or not). The order of the probability values follows the order
    the classes are given in(i.e in our example the first value is not hired and the second is hired)"""
    prob = clf.predict_proba(combinations)
    
    varprob = {}
    count = 0
    
    """Creates a dictionary of each combination and its corresponding probability"""
    for x in combinations:
        varprob[str(x)] = prob[count]
        count = count + 1
    
    return varprob