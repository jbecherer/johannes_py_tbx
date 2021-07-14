""" This package contains my implemention of an harminic analysis """


import numpy as np
import matplotlib.pylab as plt
import sys
sys.path.append('/home/johannes/Dropbox/python_tbx/johannes_py_tbx/')
import jo_plot as jpl


def my_harmo( x, y_in, F, Amin=0, Amax=2, dA=.01, Pmin=0, Pmax=2*np.pi, dP=.1):
    """ This function makes an harmonic analysis of y_in on frequencies F. 
    and returns A (amplitudes) and P (phases)"""

    # remove nans
    x = x[~np.isnan(y_in)]
    y_in = y_in[~np.isnan(y_in)] 

    # zeros element (mean)
    A_fit = np.zeros(F.shape)
    P_fit = np.zeros(F.shape)
    A_fit[0] = np.mean(y_in)
    P_fit[0] = 0

    # fitting range
    A_range = np.arange( Amin, Amax, dA)
    P_range = np.arange( Pmin, Pmax, dP)


    # remove mean
    y1 = y_in - A_fit[0]

    # error matrix

    for f in range(1,F.size):
        E = np.zeros([A_range.size, P_range.size])
        for a in range(A_range.size):
            for p in range(P_range.size):
                ytmp = A_range[a]*np.cos( - F[f]*x + P_range[p] )
                E[a,p] = np.mean( (y1-ytmp)**2 )

        min_index = np.unravel_index(E.argmin(), E.shape)
        A_fit[f] = A_range[min_index[0]]
        P_fit[f] = P_range[min_index[1]]
        y1 = y1 - A_fit[f]*np.cos( - F[f]*x + P_fit[f]  )


    return A_fit,P_fit

def get_rotaries_from_2d_field( t, U_in, F, Amin=0, Amax=2, dA=.01, Pmin=0, Pmax=2*np.pi, dP=.1):
    """ this function takes a complex velocity 2d field defined by t, U_in and
        1) makes an harmonic analysis at frequencies F in amplitude range (Amin=0, Amax=2, dA=.01)
            and Phase range (Pmin=0, Pmax=2*np.pi, dP=.1)
        2) transforms and returns the harmonics into rotay components Wac, Wc, Phi_ac, Phi_c
        """
      
    # find time dimension
    dim = np.asarray([0,1])[ np.asarray(U_in.shape)==t.size ]
    if dim[0] == 0: # if time is first dimension -> transpose
        U_in = U_in.T

    dimz = U_in.shape[0]
    dimF = F.size
    W_ac    = np.zeros((dimz, dimF))
    W_c     = np.zeros((dimz, dimF))
    Phi_ac  = np.zeros((dimz, dimF))
    Phi_c   = np.zeros((dimz, dimF))

    for i in range(dimz):
        W_ac[i,:], W_c[i,:], Phi_ac[i,:], Phi_c[i,:] =  get_rotaries_from_time_series( 
                                t, U_in[i,:], F, Amin=0, Amax=2, dA=.01, Pmin=0, Pmax=2*np.pi, dP=.1)


    return W_ac, W_c, Phi_ac, Phi_c 




def get_rotaries_from_time_series( t, U_in, F, Amin=0, Amax=2, dA=.01, Pmin=0, Pmax=2*np.pi, dP=.1):
    """ this function takes a complex velocity time series defined by t, U_in and
        1) makes an harmonic analysis at frequencies F in amplitude range (Amin=0, Amax=2, dA=.01)
            and Phase range (Pmin=0, Pmax=2*np.pi, dP=.1)
        2) transforms and returns the harmonics into rotay components Wac, Wc, Phi_ac, Phi_c
        """


    # make harmonic anaysis
    A_u, Phi_u = my_harmo( t, np.real(U_in), F, Amin, Amax, dA, Pmin, Pmax, dP)
    A_v, Phi_v = my_harmo( t, np.imag(U_in), F, Amin, Amax, dA, Pmin, Pmax, dP)


    W_ac    = np.zeros(F.shape)
    W_c     = np.zeros(F.shape)
    Phi_ac  = np.zeros(F.shape)
    Phi_c   = np.zeros(F.shape)
    for i in range(F.size):
         W_ac[i], W_c[i], Phi_ac[i], Phi_c[i] =  my_rotary_components( A_u[i], A_v[i], Phi_u[i], Phi_v[i] )


    return W_ac, W_c, Phi_ac, Phi_c 





def U_from_rotary( t, f, Wac, Wc, Phi_ac, Phi_c):
    """ This function calculates from a single rotary compoenent 
        defined by  f, Wac, Wc, Phi_ac, Phi_c
         a complex velocity vector at times t """
    U =  Wac*np.exp( 1j*f*t + 1j*Phi_ac ) + Wc*np.exp( -1j*f*t + 1j*Phi_c )
    return U


def U_from_rotaries( t, f, Wac, Wc, Phi_ac, Phi_c):
    """ This function superimposes several rotary compoenents 
        defined by  f, Wac, Wc, Phi_ac, Phi_c
        to a complex velocity vector at times t """

    U = np.zeros( t.shape)
    for i in range(f.size):
        U = U + U_from_rotary(t, f[i], Wac[i], Wc[i], Phi_ac[i], Phi_c[i])

    return U



def my_rotary_components( u0, v0, Phi_u, Phi_v ):
    """ This function returns the clock wise and anticlockwise amplitude Wac, Wc
        and phases Phi_ac, Phi_c, which together generate a complex velocity 
        vector with formular:
        U = Wac*np.exp( 1j*f*t + 1j*Phi_ac ) + Wc*np.exp( -1j*f*t + 1j*Phi_c )"""

    # van Haren 2000 version
    utild = u0*np.exp(1j*Phi_u)
    vtild = v0*np.exp(1j*Phi_v)

    wp = 0.5*(utild - 1j*vtild)
    wm = 0.5*(utild + 1j*vtild)

    Wac     =  np.abs(wp)
    Phi_ac  = -np.angle(wp)
    Wc      =  np.abs(wm)
    Phi_c   =  np.angle(wm)


    return Wac, Wc, Phi_ac, Phi_c


def my_superpos( x, f, A, P):
    """ genreates a supperposition of A[i]*cos( - f[i]*x + P[i] ) functions""" 
    y = np.zeros( x.shape)
    for i in range(f.size):
        y = y +  A[i]*np.cos( -f[i]*x + P[i] )
    return y


def coriolis_freq(lat):
    """ return f for a given latitude (lat in [deg])"""
    Omega = 7.2921e-5
    f = 2*Omega*np.sin(lat/180*np.pi)
    return f


def get_tidal_components():
    """ get a list of the main tidal compnents"""
    P = [12.42, 23.93, 12.00, 25.82, 6.21, 4.14, 24.07, 12.66, 11.97, 26.87 ]
    Labels = ['M2', 'K1', 'S2', 'O1', 'M4', 'M6',  'P1', 'N2', 'K2', 'Q1']
    Names = ['principal lunar', 'luni-solar diurnal', 'principal solar',
            'principal lunar diurnal', 'M4', 'M6' , 'principal solar diurnal', 
            'larger lunar elliptic', 'luni-solar semidiurnal', 'larger lunar elliptic']
    P = np.asarray(P)*3600
    F = 2*np.pi/P

    return F, P, Labels, Names


def test_my_harmo():
    """ example test for harmonic analysis """
    N=30
    x = np.arange( 0, N*2*np.pi, .1, dtype='f4')

    P = 2*np.pi*np.array([0, 1/4, 1/3, 0])
    f = np.array( [0, 1, 1.03, .4] )
    A = np.array( [.1, 1, .56, .3] )

    y = my_superpos( x, f, A, P ) #  + np.random.random(x.shape)

    f_fit = f[:]
    A_fit, P_fit = my_harmo( x, y, f_fit)
    y_fit = my_superpos( x, f_fit, A_fit, P_fit)

    print(A)
    print(A_fit)
    print(P)
    print(P_fit)

    plt.close('all')
    fig = plt.figure( figsize = (8, 6), facecolor = (1, 1, 1))
    ax =  jpl.create_axes(fig, 1, 1, put_xlab=False, put_ylab=False, linkx=False, linky=False)
    jpl.shift_axes(ax, 0, 0)
    jpl.squeeze_axes(ax, 1, 1)

    a=0
    ax[a].plot(x, y, label='org')
    ax[a].plot(x, y_fit, label='fit')
    ax[a].legend( loc='upper right')

    fig.show()

