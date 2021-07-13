import pandas as pd
import math
from itertools import *
import timeit
from mpmath import *
import numpy as np
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, Bounds
from scipy.special import rel_entr
import warnings
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.font_manager as font_manager

value_list = [['Male','Caucasian'],['Female','Caucasian'],['Male','African-American'],['Female','African-American'],
              ['Male','Hispanic'],['Female','Hispanic'],['Male','Asian'],['Female','Asian'],
              ['Male','Other'],['Female','Other']]

hypothesis = 10*[1/10] #uniform hypothesis

alpha = 0.05 
data = observation
bins = 1
num_bins = []
avg_run_times = []
int_avg_run_times = []

while bins < 4:
    selectedvalue = []
    for idx in range(bins):
        selectedvalue += [value_list[idx]]
    
    ndata, ncounts, nhypothesis, nvaluelist = data_slicer_multival(data, hypothesis, value_list, selectedvalue)
    my_time = %timeit -n 1 -r 10 -o hypothesis_test(ndata, nvaluelist, alpha, nhypothesis)
    #print((nvaluelist, alpha, nhypothesis))
    bins += 1
    num_bins.append(bins)
    avg_run_times.append(str(my_time))
    int_avg_run_times.append(my_time.average)

bin_count = list(range(2, 5))

#print(bin_count)
#print(avg_run_times)
run_time_frame = pd.DataFrame([bin_count, avg_run_times])
run_time_frame.to_excel(r"Time scaling experiment results.xlsx")
#graph_times(bin_count,int_avg_run_times)