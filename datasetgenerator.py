from numpy.random import default_rng
rng = default_rng()
variablelist = ['m', 'f', 'nb']
probabilitylist = [0.7, 0.2, 0.1] #0.7 corresponds to m, 0.2 to f, etc
size = 20
dataset = rng.choice(variablelist, size = size, p = probabilitylist) 
'''if you want to create paired lists (like gender+race), make a dataset for each 
individually (same size), then use zip(set1, set2) to combine them index by index'''