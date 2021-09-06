import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm
import datetime as datetime
import gsw               # sea water tool box
import jo_tools

def create_axes(fig, n_rows, n_colu, put_xlab=False, put_ylab=False, linkx=False, linky=False): 
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
            if len(ax) > 0 : # if shared axes is desired

               if linky & linkx : 
                   ax.append( fig.add_axes([(Fx[0] + (j)*(dx+xw)), (Fy[1] - (i+1)*yw - i*dy ), xw, yw],
                                           sharex=ax[0], sharey=ax[0]));
               elif linky :
                   ax.append( fig.add_axes([(Fx[0] + (j)*(dx+xw)), (Fy[1] - (i+1)*yw - i*dy ), xw, yw],
                                           sharey=ax[0]))
               elif linkx :
                   ax.append( fig.add_axes([(Fx[0] + (j)*(dx+xw)), (Fy[1] - (i+1)*yw - i*dy ), xw, yw],
                                           sharex=ax[0]))  
               else:
                   ax.append( fig.add_axes([(Fx[0] + (j)*(dx+xw)), (Fy[1] - (i+1)*yw - i*dy ), xw, yw]))

            else:
                ax.append( fig.add_axes([(Fx[0] + (j)*(dx+xw)), (Fy[1] - (i+1)*yw - i*dy ), xw, yw]))  

    
            # remove redundant lables
            if not i==(n_rows-1) and not put_xlab:
                ax[-1].tick_params(labelbottom=False)
            if not j==0 and not put_ylab:
                ax[-1].tick_params(labelleft=False)
            
            
    return ax

# }}}


def move_legend(ax, lg, dx, dy):
    # Get the bounding box of the original legend
    bb = lg.get_bbox_to_anchor().inverse_transformed(ax.transAxes)

    bb.x0 += dx
    bb.x1 += dx
    bb.y0 += dy
    bb.y1 += dy
    # Change to location of the legend.
    lg.set_bbox_to_anchor(bb, transform = ax.transAxes)


def density_contours4TS(ax):
    """This function puts lines of constant density 
    into a given T-S-asxis"""

    Tl = ax.get_ylim()
    Sl = ax.get_xlim()

    DS = (Sl[1]-Sl[0])/10
    DT = (Tl[1]-Tl[0])/10 
    s = np.arange(Sl[0], Sl[1]+DS, DS)
    t = np.arange(Tl[0], Tl[1]+DT, DT)

    S, T = np.meshgrid( s, t)
    R = gsw.sigma0(S, T)

    minR = np.min(R)
    maxR = np.max(R)
    DR = maxR-minR
    
    Rlevels = np.arange( minR, maxR, DR/7)

    conts = ax.contour(S, T, R, levels=Rlevels, colors='.7', Linewidth=1)
    #ax.clabel(conts, inline=1, fontsize=10, fmt='%3.0f m') 



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
        t = ax.text(-.01, 1.00, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner in [-2, 10]:
        t = ax.text(.5, 1.01, txt, transform=ax.transAxes, ha='center', va= 'bottom')
    elif corner == -3:
        t = ax.text(1.01, 1.0, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner == -4:
        t = ax.text(-.01, .5, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner == -5:
        t = ax.text(.5, .5, txt, transform=ax.transAxes, ha='center', va= 'center')
    elif corner == -6:
        t = ax.text(1.01, .5, txt, transform=ax.transAxes, ha='left', va= 'center')
    elif corner == -7:
        t = ax.text(-.01, -.01, txt, transform=ax.transAxes, ha='right', va= 'center')
    elif corner == -8:
        t = ax.text(.5, -.01, txt, transform=ax.transAxes, ha='center', va= 'top')
    elif corner == -9:
        t = ax.text(1.01, -.01, txt, transform=ax.transAxes, ha='left', va= 'center')
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


def reduce_panelcount(Nx, Ny, dx, dy, pos1):
  """helper function for squeese axis"""
  # double check if sum of all gaps is less than one panel width for many plots
  #import ipdb; ipdb.set_trace() # n:next line, s:step into function, c:continue, l:display position
  if (Nx-1)*dx > pos1.width:
    Nx -= 1 
    dx = (Dx-pos1.width*Nx)/(Nx-1)
  if (Ny-1)*dy > pos1.height:
    Ny -= 1 
    dy = (Dy-pos1.height*Ny)/(Ny-1)

  return Nx, Ny


def squeeze_axes(ax, Sx, Sy):
    """docstring for squeeze_axes
        this function squeezes the given axes 
        by factors Sx and Sy
    """
    if type(ax) == list:  # if more than one axes
        # fix aps between plots 
        pos1 = ax[0].get_position()
        pos2 = ax[-1].get_position()
        Dx = pos2.x0+pos2.width - pos1.x0  # total x dim of axis
        Dy = pos1.y0+pos1.height - pos2.y0  # total y dim of axis

        # number of panels in each dimension 
        # assuming all plots have the same width
        Nx = np.round(Dx/pos1.width)
        Ny = np.round(Dy/pos1.height)

        # gap size between plots
        if Nx > 1 :
          dx = (Dx-pos1.width*Nx)/(Nx-1)
        else:
          dx = 0
        if Ny > 1:
          dy = (Dy-pos1.height*Ny)/(Ny-1)
        else:
          dy = 0
        

        cnt = 1
        while Nx*Ny != len(ax):
          Nx,Ny = reduce_panelcount(Nx, Ny, dx, dy, pos1)
          cnt+=1
          if cnt>10:  # catch infinit loop
            break


        #import ipdb; ipdb.set_trace() # n:next line, s:step into function, c:continue, l:display position 
        for i in range(len(ax)):
            pos = ax[i].get_position()
            # call necessary off set
            ngapx = np.round(Nx*(pos.x0+pos.width - pos1.x0)/Dx)-1
            ngapy = np.round(Ny*(pos.y0+pos.height-pos2.y0)/Dy)-1

            ax[i].set_position([pos.x0-(ngapx*pos.width*(1-Sx)), pos.y0-(ngapy*pos.height*(1-Sy)), pos.width * Sx , pos.height * Sy])

    else:
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width * Sx , pos.height * Sy])
        

def copy_axes():
    """docstring for copy_axes"""
    pass

def tstamp2yday(t):
    """ converts time stamp into yday """
    yday = jo_tools.tstamp2yday( t )
    return yday

def mypcolor(ax, x, y, z, cmap=cm.RdYlBu, cl=None, colbar=True, title=None, ylab=None):
    """ quick pcolor function for a given axis including colorbar and title"""

    if x[0] > 1e9:   # if x is a unix time stemp convert into year days
        x = tstamp2yday(x)

    if cl:
        cl = np.asarray(cl)
        pc = ax.pcolor(x, y, z, cmap=cmap , vmin=cl[0], vmax=cl[1])	
    else:
        pc = ax.pcolor(x, y, z, cmap=cmap) # , vmin=0, vmax=1e-2)	

    if colbar:
        fig = ax.get_figure()
        axpos = ax.get_position()
        cax  = fig.add_axes([axpos.x1+.02, axpos.y0+.1*axpos.height, .01, .8*axpos.height])
        cb = fig.colorbar(pc , extend='both', cax=cax) # ticks=[1,2,3])
    else:
        cb = []

    if title:
        txt = text_corner(ax, title, 1)
        txt.set_backgroundcolor([1, 1, 1, .5])
        txt.set_fontsize(12)
        txt.set_color([0, 0, 0])
    else:
        txt = ''

    if ylab:
        ax.set_ylabel(ylab, fontsize=12)

    return pc, cb, txt


def plot_x_line(ax, val, style='k--') :
  """ plot horizontal line in ax at value val with style"""

  xl = ax.get_xlim()
  ax.plot(xl, np.array([1,1])*val, style)
  ax.set_xlim(xl)
  


# Test create_axes()
if False :
    fig = plt.figure( figsize = (10, 10), facecolor = (1, 1, 1))
    ax =  create_axes(fig, 3, 2, put_xlab=False, put_ylab=False)
    shift_axes(ax[0], -.03, .03)
    shift_axes(ax[slice(2,-1,2)], -.01, -.03)
    squeeze_axes(ax[slice(2,-1,2)], .8, 1.05)
    ax[-1].text( 1, 1, str(1), va='top', ha='right')
    t = text_corner(ax[-1], 'corner text', 2)
    t = text_corner(ax[-1], 'corner text', 1245)

    abc = 'abcdefghijklmnopqrstuvwxyz'
    for a in range(len(ax)):
        tabc = text_corner(ax[a], '(' + abc[a] + ')', 7);
        tabc.set_fontweight('bold')
        tabc.set_backgroundcolor([1, 1, 1, .5])

    fig.show() 

