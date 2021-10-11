import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm
import os
import sys
import datetime as datetime
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

def movmedian(X, n=1, axis=0, maskNAN=False):
    """This function generates a moving median ignoring nan along the axis dimension 
    of X of oder n"""


    ndim = X.ndim

    if maskNAN:
        mnan = np.isnan(X)

    if n==0:
        return X

    if axis==1: # if flipped dimensions
        Xout = np.transpose( movmedian( np.transpose(X), n))
        return Xout

    if n>1: # itteration for higher order average
        X = movmedian(X, n=n-1)

    if ndim == 2 :
        xleft  = np.full( (X.shape[0], X.shape[1]), np.nan, dtype='f4')
        xright = np.full( (X.shape[0], X.shape[1]), np.nan, dtype='f4')
        xleft[:-1,:] = X[1:,:] 
        xright[1:,:] = X[:-1,:]
        Xout = np.nanmedian( np.concatenate( 
                        (np.tile(X,(1,1,1)), np.tile(xleft,(1,1,1)), np.tile(xright,(1,1,1))),
                            axis=0), axis=0)
    else:
        xleft  = np.full( (X.shape[0], ), np.nan, dtype='f4')
        xright = np.full( (X.shape[0], ), np.nan, dtype='f4')
        xleft[:-1] = X[1:] 
        xright[1:] = X[:-1]
        Xout = np.nanmedian( np.concatenate( 
                        (np.tile(X,(1,1)), np.tile(xleft,(1,1)), np.tile(xright,(1,1))),
                            axis=0), axis=0)

    if maskNAN:
        Xout[mnan] = np.nan

    return Xout


def interpNans(x):
  """This function interps over nans in signal 
    OUTPUT: x_new
  """
  ii_nan = np.where(np.isnan(x))[0]
  ii_notnan = np.where(~np.isnan(x))[0]
  x_new = x.copy()
  x_new[ii_nan] = np.interp( ii_nan, ii_notnan, x[ii_notnan])

  return x_new

def interpWithNans(xnew, xold, yold):
  """This function interps despite nans in signal 
    OUTPUT: ynew
  """
  ii_notnan = np.where( ~(np.isnan(xold) | np.isnan(yold)) )[0]

  ynew = np.full( xnew.shape, np.nan, dtype='f8')

  if ii_notnan.size > 1 :
    ynew = np.interp( xnew, xold[ii_notnan],  yold[ii_notnan])

  return ynew


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


def ll_dist( lat1, lon1, lat2, lon2 ):
  """ calculates distance (km) between two points based on the 
  assumption of a spherical Earth
  """
  R = 6373.0

  lat1 = np.radians(lat1)
  lon1 = np.radians(lon1)
  lat2 = np.radians(lat2)
  lon2 = np.radians(lon2)

  dlon = lon2 - lon1
  dlat = lat2 - lat1
  a = (np.sin(dlat/2))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2))**2
  c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
  distance = R * c

  return distance


def tstamp2yday( tstep ):
  """ This function converts unix time to year day"""
  # day of the year
  day1year = datetime.datetime(year=datetime.datetime.utcfromtimestamp(tstep[1]).year, month=1, day=1, 
      tzinfo=datetime.timezone.utc)
  yday = (tstep - day1year.timestamp())/3600/24 + 1
  
  return yday

def yday2tstamp( yday, year ):
  """ converts yday to unix time """
  day1year = datetime.datetime(year=year, month=1, day=1, tzinfo=datetime.timezone.utc)
  tstamp = (yday-1)*3600*24 + day1year.timestamp()
  
  return tstamp

def tstamp2dtime( tstamp ):
  """ converts unix time to date time object """
  
  dtime = [datetime.datetime.utcfromtimestamp(t) for t in tstamp]
  
  return dtime

def yday2dtime( yday, year ):
  """ converts year day to date time object """
  tstamp = yday2tstamp( yday, year )
  dtime  = tstamp2dtime( tstamp )
  return dtime


def bin_interp( x, y, xnew):
  """This function first binaverages x,y in xnew
     and then interpolates over empty bins"""

  # construct bins with xnew as center points
  bins = x2bins(xnew)

  ybin_avg = binavg( x, y, bins)

  # interp over nans in bins
  ii_nan = np.isnan(ybin_avg)
  ii_notnan = ~np.isnan(ybin_avg)
  ynew = ybin_avg.copy()
  if  (np.where( ii_nan )[0].size > 0) & (np.where( ii_notnan )[0].size > 0) : 
    ynew[ii_nan] = np.interp( xnew[ii_nan], xnew[ii_notnan], ynew[ii_notnan])

  return ynew



def binavg( x, y, bins):
  """This function first binaverages x,y in bins
     """

  inds = np.digitize(x, bins)
  ybin_avg = np.array([ y[inds == i].mean() for i in range(1, len(bins))])

  return ybin_avg
  
def x2bins(x):
  """This function construct bins with x as center point"""

  bins = np.zeros( (x.size+1, ) )
  bins[0] = x[0] - .5*(x[1]-x[0])
  bins[-1] = x[-1] + .5*(x[-1]-x[-2])
  bins[1:-1] = .5*(x[0:-1] + x[1:])

  return bins

