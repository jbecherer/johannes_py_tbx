import numpy as np
import matplotlib.pylab as plt

def create_axes(fig, n_rows, n_colu, put_xlab=False, put_ylab=False): 
    """ # {{{
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

# }}}

def text_corner(ax, txt, corner): 
    """ # {{{ 
        t = text_corner(ax, txt, corner)  
         this function generates a text (txt) in the corner of the 
         axis (ax) , where (t) is the text property handel
         t = text_corner(ax, txt, corner)"""
    if corner == 1:
        t = ax.text(.01, .99, txt,  transform=ax.transAxes,ha='left', va= 'top')
    elif corner == 2:
        t = ax.text(.5, .99, txt, transform=ax.transAxes, ha='center', va='top')
    elif corner == 3:
        t = ax.text(.99, .99, txt, transform=ax.transAxes, ha='right', va='top',)
    elif corner == 4:
        t = ax.text(.01, .5, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner == 5:
        t = ax.text(.5, .5, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner == 6:
        t = ax.text(.99, .5, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner == 7:
        t = ax.text(.01, .01, txt, transform=ax.transAxes, ha='left', va= 'bottom')
    elif corner == 8:
        t = ax.text(.5, .01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner == 9:
        t = ax.text(.99, .01, txt, transform=ax.transAxes, ha='right', va= 'bottom')
    elif corner == 0:
        t = ax.text(.5, 1.01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner == -1:
        t = ax.text(.01, 1.01, txt, transform=ax.transAxes, ha='left', va= 'bottom')
    elif corner in [-2, 10]:
        t = ax.text(.5, 1.01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner == -3:
        t = ax.text(.99, 1.01, txt, transform=ax.transAxes, ha='left', va= 'bottom')
    elif corner == -4:
        t = ax.text(-.01, .5, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner == -5:
        t = ax.text(.5, .5, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner == -6:
        t = ax.text(1.01, .5, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner == -7:
        t = ax.text(.01, -.01, txt, transform=ax.transAxes, ha='left', va= 'top')
    elif corner == -8:
        t = ax.text(.5, -.01, txt, transform=ax.transAxes, ha='center', va= 'top')
    elif corner == -9:
        t = ax.text(.99, -.01, txt, transform=ax.transAxes, ha='right', va= 'top')
    elif corner in [12 , 21]:
        t = ax.text(.25, .99, txt, transform=ax.transAxes, ha='center', va= 'top')
    elif corner in [23 , 32]:
        t = ax.text(.75, .99, txt, transform=ax.transAxes, ha='center', va= 'top')
    elif corner in [45 , 54]:
        t = ax.text(.25, .5, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [56, 65]:
        t = ax.text(.75, .5, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [78 , 87]:
        t = ax.text(.25, .01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner in [89 , 98]:
        t = ax.text(.75, .01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner in [14, 41]:
        t = ax.text(.01, .75, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner in [47, 74]:
        t = ax.text(.01, .25, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner in [36, 63]:
        t = ax.text(.99, .75, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner in [69, 96]:
        t = ax.text(.99, .25, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner in [25, 52]:
        t = ax.text(.5, .75, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [58, 85]:
        t = ax.text(.5, .25, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [1245, 4512, 1425, 2514]:
        t = ax.text(.25, .75, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [2356, 5623, 2536, 3625]:
        t = ax.text(.75, .75, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [4578, 7845, 4758, 5847]:
        t = ax.text(.25, .25, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [5869, 5689, 6958, 8956, 9865]:
        t = ax.text(.75, .25, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner in [-12 , -21]:
        t = ax.text(.25, 1.01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner in [-23 , -32]:
        t = ax.text(.75, 1.01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner in [-78 , -87]:
        t = ax.text(.25,-.01, txt, transform=ax.transAxes, ha='center', va= 'top')
    elif corner in [-89 , -98]:
        t = ax.text(.75,-.01, txt, transform=ax.transAxes, ha='center', va= 'top')
    elif corner in [-14, -41]:
        t = ax.text(-.01, .75, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner in [-47, -74]:
        t = ax.text(-.01, .25, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner in [-36, -63]:
        t = ax.text(1.01, .75, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner in [-69, -96]:
        t = ax.text(1.01, .25, txt, transform=ax.transAxes, ha='left', va= 'center')
    else:
        print('choose corner:');
        print('type text_corner() without arguments to see options');
        return []
    
    return t   # }}}

def shift_axes(ax, dx, dy):
    """docstring for shift_axes%            
        this function shifts a given set of axes ax 
        by dx and dy """
    if type(ax) == list:  # if more than one axes
        for a in range(len(ax)):
            pos = ax[a].get_position()
            ax[a].set_position([pos.x0 + dx, pos.y0 + dy, pos.width, pos.height])
    else:
        pos = ax.get_position()
        ax.set_position([pos.x0 + dx, pos.y0 + dy, pos.width, pos.height])



def squeeze_axes():
    """docstring for squeeze_axes"""
    pass

def copy_axes():
    """docstring for copy_axes"""
    pass



# Test create_axes()
fig = plt.figure( figsize = (10, 10), facecolor = (1, 1, 1))
ax =  create_axes(fig, 3, 2, put_xlab=False, put_ylab=False)
shift_axes(ax[0], -.03, .03)
shift_axes(ax[slice(2,-1,2)], -.01, -.03)
ax[-1].text( 1, 1, str(1), va='top', ha='right')
t = text_corner(ax[-1], 'corner text', 2)
t = text_corner(ax[-1], 'corner text', 1245)

abc = 'abcdefghijklmnopqrstuvwxyz'
for a in range(len(ax)):
    tabc = text_corner(ax[a], '(' + abc[a] + ')', 7);
    tabc.set_fontweight('bold')
    tabc.set_backgroundcolor([1, 1, 1, .5])



fig.show() 

