import numpy as np
from Initialize import *
print(F21)

## Defining functions, 
def emissions(yr, e_stop = True):
    """Function defining A2 emission scenario over the interval 1990-2100
    extended to pre-industrial (assuming linear increase from 0 in 1850 to
    1900) and assuming full cessation of CO_2 input at 2101
    For additional information see http://www.grida.no/climate/ipcc/emission

    Parameters
    ----------
    yr: 1D array
        year range for emission scenario

    Returns
    -------
    e: double
        A2 emissions in a given year

    """

    t_yr = np.array([0, 1850, 1990, 2000, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120, 10000])
    e_GtC_yr = np.array([0, 0,  6.875, 8.125, 9.375, 12.5, 14.375, 16.25, 17.5, 19.75, 21.25, 23.125, 26.25, 28.75, 0, 0, 0])

    if e_stop == False:
        # changed last 3 variables in vector to the 2000 year emission value (previously all 0)
        e_GtC_yr = np.array([0, 0,  6.875, 8.125, 9.375, 12.5, 14.375, 16.25, 17.5, 19.75, 21.25, 23.125, 26.25, 28.75, 8.125, 8.125, 8.125])

    e = np.interp(yr, t_yr, e_GtC_yr)

    return e

def Find_k(F_in, F_out, M):
    """
    Finds values of linear flux constants, k, for fluxes in and out for a given set of fluxes in, out and inital masses. 
    ## Inputs ##
    F_in: NxM array
        Array with flux-in values
    F_out: NxM array
        Array with flux-out values, positive flux values
    M: M length vector
        Initial values of mass for each box
    ## Outputs ##
    k: NxM array
        Array of k constant values to solve M-box model. k*M gives M ODEs to solve M box model problem
    """
    k_in = F_in.T/M
    k_out = sum(F_out/M)
    k = k_in.T - k_out*np.diag(np.ones(len(F_in)))
    return k

def MassFlux_forced(t, M):
    """
    Mass flux coupled ODEs for 4-box model of carbon cycle in a steady state
    """
    ## Flux in and Flux out values, units Gt/yr, to find k values
    if len(M) == 4:
        k = Find_k(Flux_in_4, Flux_out_4, M0_4)
    else:
        k = Find_k(Flux_in_9, Flux_out_9, M0_9)

    ## Matrix multiplication to create linear system of equations
    dMdt = np.matmul(k, M) # this result is the steady state solution
    dMdt[0] += emissions(t) ## adding forcing from emissions function at each timestep
    return dMdt

def MassFlux_Steady(t, M):
    """
    Unforced coupled ODEs for 4-box model of carbon cycle. 
    """
    global k_4, k_9
    ## Flux in and Flux out values, units Gt/yr, to find k values
    if len(M) == 4:
        k = Find_k(Flux_in_4, Flux_out_4, M0_4)
    else:
        k = Find_k(Flux_in_9, Flux_out_9, M0_9)

    ## Matrix multiplication to create linear system of equations
    dMdt = np.matmul(k, M) # this result is the steady state solution
    return dMdt

def rk4(fxy, x0, xf, y0, N):
    ####### From Mark (or luke? Not my work as you can tell) #######
    # The inputs to the function are:
    #         fxy = the name of the function containing f(x,y) (e.g. oneode, twoode)
    #         xo,xf = initial and final values of the independent variable (integers or floats)
    #         yo = initial value of dependent variable at xo (numpy array)
    #         N = number of intervals to use between xo and xf (integer)

    # The outputs to the function are:
    #         X = numpy array containing values of the independent variable
    #         Y = the estimated dependent variable at each value of the independent variable
    #         --> this variable is a 1D numpy array if only one equation is solved
    #         --> it is an M-D numpy array [y1(x) y2(x) ... ] for multiple (M) equations 

    #compute step size and size of output variables
    if N < 2:
        N = 2 #set minimum number for N
    h = (xf - x0) / N
    X = np.zeros((N+1, 1))
    P = np.max(np.shape(y0))
    Y = np.zeros((N+1, P))*1j #make complex by multiplying by 1j; this way can add complex values to this during integration

    #set initial conditions
    x = x0
    X[0] = x
    y = [complex(val) for val in y0]  #make complex
    Y[0,:] = y
    
    #begin computational loop
    for ii in range(N):

        k1 = np.array([h * val for val in fxy(x,y)]) #evaluate function fxy; depending on equation, k1-4 can be complex; this is why we make Y and y complex as well
        k2 = np.array([h * val for val in fxy(x+h/2, y+k1/2)])
        k3 = np.array([h * val for val in fxy(x+h/2, y+k2/2)])
        k4 = np.array([h * val for val in fxy(x+h, y+k3)])

        y += (k1 + 2*k2 + 2*k3 + k4) / 6.
        x += h
        X[ii+1] = x
        Y[ii+1,:] = y
    
    return X, Y
