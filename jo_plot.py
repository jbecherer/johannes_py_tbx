import numpy as np
import matplotlib.pylab as plt

def create_axes(fig, n_rows, n_colu, put_xlab=False, put_ylab=False): # {{{
    """
    ax = create_axes(fig, n_rows, n_colums, put_xlab = False, put_ylab = False)
        this function generates for a given figure "fig" 
        a matrix of subplots "ax" with n_rows and n_colums

        2020-06-24
        Johannes Becherer
    """
    
    dx  = .015  # lateral distance between plots
    dy  = .015  # vertical sitance between plots
    
    # if labels are between panels
    if put_xlab:
        dy = 0.08
    if put_ylab:
        dx = 0.08

    # generate  frame inside the figure
    Fx = [.1, .95] # x-dimension of frame
    Fy = [.1, .95] # y-dimension of frame   



    xw = ((Fx[1] - Fx[0]) - (n_colu-1)*dx)/n_colu  # width of each subplot
    yw = ((Fy[1] - Fy[0]) - (n_rows-1)*dy)/n_rows      # width of each subplot


    ax = [] 
    for i in range(n_rows):
        for j in range(n_colu):
            ax.append( fig.add_axes([(Fx[0] + (j)*(dx+xw)), (Fy[1] - (i+1)*yw - i*dy ), xw, yw]));  
    
            # remove redundant lables
            if not i==(n_rows-1) and not put_xlab:
                ax[-1].tick_params(labelbottom=False)
            if not j==0 and not put_ylab:
                ax[-1].tick_params(labelleft=False)
            
            
    return ax

## Test create_axes()
#fig = plt.figure( figsize = (10, 10), facecolor = (1, 1, 1))
#ax =  create_axes(fig, 10, 10, put_xlab=False, put_ylab=False)
#ax[-1].text( 1, 1, str(1), verticalalignment='top', horizontalalignment='right')
#fig.show() 
# }}}

