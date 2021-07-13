while bins < 10:
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

bin_count = list(range(2, 11))

#print(bin_count)
#print(avg_run_times)
run_time_frame = pd.DataFrame([bin_count, avg_run_times])
run_time_frame.to_excel(r"Time scaling experiment results.xlsx")