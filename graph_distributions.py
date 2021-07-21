
def graph_distributions(q,value_list, hypothesis=[],selected_value=[],filename="distributions.pdf"):

    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib import cm
    import numpy as np
    import math
    
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
        
    nvalue_list = []

    if len(selected_value) != 0:
        for x in range(len(selected_value)):
            #flat = [val for sublist in x for val in sublist]
            nvalue_list.append(selected_value[x])
            #print(string)
            not_val = ['NOT']*1 + selected_value[x]
            nvalue_list.append(not_val)
            #nvalue_list = [[x in selected_value], ['NOT MALE WHITE']]
        print(nvalue_list)
    else:
        nvalue_list = value_list

    N = len(nvalue_list)
    #print(N)

    if (N > 11):
        text_size = 3
        h_line_factor = .2
        xy_factor = 8
    elif(N < 5 and N >= 2):
        text_size = 5
        h_line_factor = .05
        xy_factor = 11
    else:
        text_size = 5
        h_line_factor = .2
        xy_factor = 11
    
    #print(text_size,h_line_factor,xy_factor)

    xl = []
    for x in nvalue_list:
        header = ''
        num_classes = len(x)
        for i in range(num_classes):
            n = x[i]
            if n.count('-') > 0:
                n = n.replace("-"," ")
            if n.count(' ') > 0:
                split = [s[0] for s in n.split()]
                for ind in range(len(split)):
                    header = str(header) + str(split[ind])
            else:
                header = str(header) + str(str(n[:1]))
        xl.append(header)
    #print(xl)

    xmin, xmax = xlim = 0, N+1
    ymin, ymax = ylim = 0, 1.1

    fig, ax = plt.subplots(dpi=300)
    ax.set(xlim=xlim, ylim=ylim, autoscale_on=False)


    #xl = ['MW','FW','MOA','FOA','MAA','FAA','MOH','FOH','MA','FA']
    x = np.arange(1,N+1) + .15
    yh = np.array(hyp)#np.random.rand(N)
    bar_width = 0.1
    y = np.array(q)#np.random.rand(N)


    color_scheme = matplotlib.cm.get_cmap('inferno')
    norm = plt.Normalize(-1,1)
    sm = cm.ScalarMappable(cmap=color_scheme,norm=norm)  # norm sets the maximum and minimum values
    sm.set_array([])
    plt.colorbar(sm)

    for idx in range(N):

        diff = 0.5 + (y[idx]-yh[idx])
        
        
        plt.hlines(y=yh[idx],xmin=(x[idx]-h_line_factor),xmax=(x[idx]+h_line_factor), colors='black',alpha=0.3,
                       linestyle='solid',linewidth = 1, zorder=3)
       
        plt.vlines(x[idx], ymin=0, ymax=y[idx], colors=color_scheme(diff), linestyles='solid',linewidth=3)
        plt.plot(x[idx],y[idx],"o",color=color_scheme(diff))
       
        plt.annotate(format(y[idx], '.2f'), 
                                (x[idx], y[idx]), 
                                ha = 'center', va = 'center',
                                xytext = (0, 5.5), 
                                textcoords = 'offset points', weight='bold',fontsize=text_size)
       
        
             
           
       


    legend_elements = [Line2D([0], [0], color='black',alpha=0.3, linestyle='solid',lw=1, label='Proposed'),
                       Line2D([0], [0], color='w', marker='o',markerfacecolor='black',markersize=6,
                              lw=4, label='Plausible')]

    ax.set_aspect('auto')
    ax.set_xticks(x)
    ax.set_xticklabels(xl,rotation=30)
    ax.set_ylabel("Proportion",fontsize = 11,fontname="Sans-serif")
    ax.set_title("Proposed Distribution vs. Closest Plausible Distribution",fontname="Sans-serif")
    ax.legend(handles=legend_elements,fontsize=6,loc='upper right')
    plt.savefig(filename,bbox_inches="tight",pad_inches=0) #saves the image as a pdf. you can change the file name to any name of format you want. this is the default I came up with.
