import matplotlib.pyplot as plt
import matplotlib
import math
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm
import numpy as np

def graph_distributions(q, value_list, hypothesis=[], selected_value=[], xl=[], filename="distributions.pdf"):
    """
    This function is used to generate the distribution comparison graphs seen
    in the Results section of the paper.
    
    ***The function is only setup to graph up to 10 values in a value list. However, bar 
    spacing, text, and hlines can be edited to accomodate larger value list***
    """
    
    # Sets hypothesis to uniform hypothesis if one is not given
    if hypothesis == []:
        hyp = len(value_list)*[1/len(value_list)]
    else:
        hyp = hypothesis
        
    nvalue_list = []
    
    # Create New value list called nvalue_list which puts value_list in acceptable
    # format for graph based on given value_list
    if len(selected_value) != 0:
        for x in range(len(selected_value)):
            nvalue_list.append(selected_value[x])
            not_val = ['NOT']*1 + selected_value[x]
            nvalue_list.append(not_val)
        print(nvalue_list)
    else:
        nvalue_list = value_list

    N = len(nvalue_list) # N denotes the number of values being graphed
    
    # Sets text size, length of horizontal line based on number of values graphed
    if(N <= 5):
        text_size = 10
        h_line_factor = .1
    else:
        if N == 6:
            text_size = 8
        elif N == 7 or N == 8:
            text_size = 7
        elif N == 9 or N == 10:
            text_size = 6
        h_line_factor = .2
        
    # Create Xtick labels based on given nvalue_list
    # If value is hyphenated or has multiple words with spaces (i.e. 'African-American' or 'White Male'), 
    # the first letter of each word is used to create abrev. for value
    # If the nvalue_list value is a list or tuple then he first letter of each word in the list
    # is used to create abrev. for value
    if xl != []:
        for x in nvalue_list:
            header = ''
            if type(x) is list or type(x) is tuple:
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
            else:
                if x.count('-') > 0:
                    x = x.replace("-"," ")
                if x.count(' ') > 0:
                    words = x.split()
                    split = [word[0] for word in words]
                    for ind in range(len(split)):
                        header = str(header) + str(split[ind])
                else:
                    header = x
                xl.append(header)
    else:
        xl = xl
    
    # Create Ymax for graph
    if max(q) > max(hyp):
        factor = 10.0 ** 1
        yrange = math.trunc((max(q)+.2) * factor) / factor
    else:
        factor = 10.0 ** 1
        yrange = math.trunc((max(hyp)+.2) * factor) / factor
    
    # Sets Yticks to even proportion values between 0 and 1
    yticks=[0]
    tick = 0
    while tick <= yrange:
        tick = max(yticks)+.2
        
        if tick >= 1:
            break
        else:
            yticks.append(tick)

    # Set X and Y range
    xmin, xmax = xlim = 0, N+1
    ymin, ymax = ylim = 0, yrange

    fig, ax = plt.subplots(dpi=300)
    ax.set(xlim=xlim,ylim=ylim, autoscale_on=False)

    # Set xtick, y values, lollipop width, and q list
    x = np.arange(1,N+1) + .15
    yh = np.array(hyp)
    bar_width = 0.1
    y = np.array(q)

    # Create Color scheme and Colorbar
    color_scheme = matplotlib.cm.get_cmap('inferno')
    norm = plt.Normalize(-1,1)
    sm = cm.ScalarMappable(cmap=color_scheme,norm=norm)  # norm sets the maximum and minimum values
    sm.set_array([])
    plt.colorbar(sm)

    # Create Lollipop graph with q distribution threshold line for each lollipop
    for idx in range(N):
        diff = 0.5 + (y[idx]-yh[idx])
        plt.hlines(y=yh[idx],xmin=(x[idx]-h_line_factor),xmax=(x[idx]+h_line_factor), colors='black',alpha=0.3,\
                   linestyle='solid',linewidth = 1, zorder=3)
        plt.vlines(x[idx], ymin=0, ymax=y[idx], colors=color_scheme(diff), linestyles='solid',linewidth=3)
        plt.plot(x[idx],y[idx],"o",color=color_scheme(diff))
        plt.annotate(format(y[idx], '.2f'), 
                     (x[idx], y[idx]), 
                     ha = 'right', va = 'center',
                     xytext = (-5.5, 0), 
                     textcoords = 'offset points', weight='bold',fontsize=text_size)  
        
    legend_elements = [Line2D([0], [0], color='black',alpha=0.3, linestyle='solid',lw=1, label='Proposed'),
                       Line2D([0], [0], color='w', marker='o',markerfacecolor='black',markersize=6,
                              lw=4, label='Plausible')]

    ax.set_aspect('auto')
    ax.set_xticks(x)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xl,rotation=30)
    ax.set_ylabel("Proportion",fontsize = 11,fontname="Sans-serif")
    ax.set_title("Proposed vs. Closest Plausible",fontname="Sans-serif",weight='bold')
    ax.legend(handles=legend_elements,fontsize=6,loc='upper right')
    plt.savefig(filename,bbox_inches="tight",pad_inches=0) #saves the image as a pdf. you can change the file name to any name of format you want. this is the default I came up with.