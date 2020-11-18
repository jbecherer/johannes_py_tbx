import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm
import os
import sys
import datetime as dt
import scipy.io
from scipy import signal
from importlib import reload
sys.path.append('/home/johannes/Dropbox/python_tbx/johannes_py_tbx/')
import jo_plot as jpl


def gaussfilter(X, n=1):
    """makes a 2d gause filter of order n for 2d matrix X"""

    if n==0:
        return X

    X = movmean(X, n=n)
    Xout = movmean(X, n=n, axis=1)

    return Xout


def movmean(X, n=1, axis=0, maskNAN=False):
    """This function generates a moving mean ignoring nan along the axis dimension 
    of X of oder n"""


    ndim = X.ndim

    if maskNAN:
        mnan = np.isnan(X)

    if n==0:
        return X

    if axis==1: # if flipped dimensions
        Xout = np.transpose( movmean( np.transpose(X), n))
        return Xout

    if n>1: # itteration for higher order average
        X = movmean(X, n=n-1)

    if ndim == 2 :
        xleft  = np.full( (X.shape[0], X.shape[1]), np.nan, dtype='f4')
        xright = np.full( (X.shape[0], X.shape[1]), np.nan, dtype='f4')
        xleft[:-1,:] = X[1:,:] 
        xright[1:,:] = X[:-1,:]
        Xout = np.nanmean( np.concatenate( 
                        (np.tile(X,(1,1,1)), np.tile(xleft,(1,1,1)), np.tile(xright,(1,1,1))),
                            axis=0), axis=0)
    else:
        xleft  = np.full( (X.shape[0], ), np.nan, dtype='f4')
        xright = np.full( (X.shape[0], ), np.nan, dtype='f4')
        xleft[:-1] = X[1:] 
        xright[1:] = X[:-1]
        Xout = np.nanmean( np.concatenate( 
                        (np.tile(X,(1,1)), np.tile(xleft,(1,1)), np.tile(xright,(1,1))),
                            axis=0), axis=0)

    if maskNAN:
        Xout[mnan] = np.nan

    return Xout





def qbutter(x, coff, btype='low'):  # {{{
   """this function generates a filtered version of the input
   x, with cutoff coff and filter type btype='low' 
   other options are 'high'
   or 'band'.... for band you need coff=[.1,.8] for band pass
   """


   if type(coff)==list:

       if len(coff)==2:
           if coff[0]>coff[1]:
               ctmp = coff.copy()
               coff[0] = ctmp[1]
               coff[1] = ctmp[0]
        
       if coff[0]>=1.0:
            return x
   else:
       if coff>=1.0:
            return x

   if x.ndim>1:
      x = x.squeeze()


   b, a = signal.butter(3, coff, btype=btype)


   xfilt = signal.filtfilt(b, a, x)

   return xfilt # }}}


def find_clusters(time, maxdt):
    """This function finds clustered data in iregular but consequtive time series
        Usage:
            ind_ss = jo_tools.find_clusters( np.asarray(a), 1)
            [a[ ind_ss[i,0] : ind_ss[i,1] ] for i in range(ind_ss.shape[0])]
        """
    diff_time = np.diff(time)
    ind_jumps = np.where( diff_time>maxdt )

    N_jumps = ind_jumps[0].size
    ind_ss = np.full( (N_jumps+1,2), np.nan, dtype='i4')
 
    ind_ss[0,0] = 0;
    for i in range(N_jumps):
        ind_ss[i,1] = ind_jumps[0][i] + 1
        ind_ss[i+1,0] = ind_ss[i,1]
     
    # final stop
    ind_ss[-1,1] = time.size

    return ind_ss
    

def average_cluster( time, data, maxdt): # {{{
    """This function finds all clustered data and averages them together"""

    cluster_time = []
    cluster_mean = []
    cluster_median = []
    cluster_min  = []
    cluster_max  = [] 
    cluster_std  = []
    cluster_N    = []

    diff_time = np.diff(time)
    ind_jumps = np.where( diff_time>maxdt )
    
    
    N_jumps = ind_jumps[0].size
    istart = 0
    for i in range(N_jumps):
        istop = ind_jumps[0][i] + 1
        cluster_time.append( np.mean(time[istart:istop]) )
        cluster_mean.append( np.mean(data[istart:istop]) )
        cluster_median.append( np.median(data[istart:istop]) )
        cluster_min.append( np.min(data[istart:istop]) )
        cluster_max.append( np.max(data[istart:istop]) )
        cluster_std.append( np.std(data[istart:istop]) )
        cluster_N.append(istop-istart)

        istart = istop

    # last cluster
    istop = time.size
    cluster_time.append( np.mean(time[istart:istop]) )
    cluster_mean.append( np.mean(data[istart:istop]) )
    cluster_median.append( np.median(data[istart:istop]) )
    cluster_min.append( np.min(data[istart:istop]) )
    cluster_max.append( np.max(data[istart:istop]) )
    cluster_std.append( np.std(data[istart:istop]) )
    cluster_N.append(istop-istart)


    # put data to dictonary
    cluster = dict( time=cluster_time,
                     mean=cluster_mean,       
                     median=cluster_median,
                     min=cluster_min,
                     max=cluster_max,
                     std=cluster_std,
                     N=cluster_N)

    return cluster # }}}
