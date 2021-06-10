import math

def abs_sum_mg(obs_bin,hypothesis):
    """takes in both observation and hypothesis in bin form to simplify inputs. 
    For instance a uniform hypothesis would be something like [3,3,3,3].
    Interestingly, (binA, binB) =/= (binB, binA)"""
    obs_distance = [obs_bin[i]-hypothesis[i] for i in range(len(obs_bin))]
    obs_abs_sum = sum([abs(val) for val in obs_distance])
    total = sum(obs_bin)
    max_abs_sum = (total - min(hypothesis))*2
    mg = 0
    for i in range(obs_abs_sum, max_abs_sum, 2):
        pm_pair = [i//2,-i//2]
        bad_list = helper_abs_sum(hypothesis,pm_pair)
        #i don't want to have to filter but i'm not sure how to fix my helper function
        true_list = list(filter(lambda x: True if sum(x)==0 else False, bad_list))
        #print(true_list)
        for perm in true_list:
            event = [hypothesis[i]+perm[i]for i in range(len(perm))]
            mg += math.factorial(total) // math.prod(list(map(lambda x: math.factorial(x), event)))
    return mg

def helper_abs_sum(hyp,pm_pair):
    if len(hyp) == 1:
        if pm_pair[0] == 0:
            return [[max(pm_pair[1],-1*hyp[0])]]
        else:
            return [[pm_pair[0]]]
    else:
        permlist = []
        vallist = list(range(max(-1*hyp[0],pm_pair[1]), pm_pair[0]+1))
        for val in vallist:
            if abs(val) == val:
                permlist += [[val]+perm for perm in helper_abs_sum(hyp[1:], [pm_pair[0]-val, pm_pair[1]])]
            else:
                permlist += [[val]+perm for perm in helper_abs_sum(hyp[1:], [pm_pair[0], pm_pair[1]-val])]
        return permlist



            
