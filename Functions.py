import numpy as np
import matplotlib.pyplot as plt
from Initialize import *

## Defining functions, 
def emissions(yr, opt = 100):
    """Function defining A2 emission scenario over the interval 1990-2100
    extended to pre-industrial (assuming linear increase from 0 in 1850 to
    1900) and assuming full cessation of CO_2 input at 2101
    For additional information see http://www.grida.no/climate/ipcc/emission

    Parameters
    ----------
    yr: 1D array
        Year range for emission scenario

    Returns
    -------
    e: double
        A2 emissions in a given year

    """
    t_yr = np.array([0, 1850, 1990, 2000, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110, 2120, 10000])
    e_GtC_yr = np.array([0, 0,  6.875, 8.125, 9.375, 12.5, 14.375, 16.25, 17.5, 19.75, 21.25, 23.125, 26.25, 28.75, 0, 0, 0])

    period = 70 #in years
    amp = 10 # amplitude of periodic function
    if opt == 0: # changed last 3 variables in vector to the 2000 year emission value (previously all 0)
        e_GtC_yr = np.array([0, 0,  6.875, 8.125, 9.375, 12.5, 14.375, 16.25, 17.5, 19.75, 21.25, 23.125, 26.25, 28.75, 8.125, 8.125, 8.125])
    elif opt == 1: # periodic emissions (sin)
        e_GtC_yr = amp*np.sin(((2*np.pi)/period)*(t_yr - 1850)) + e_GtC_yr
    elif opt == 2: # periodic decay of emissions
        e_GtC_yr = amp*np.exp(-((2*np.pi/(10*period))*(t_yr - 2100))) + e_GtC_yr # period multiplication 10* is all exp, 100* ~same as no forcing
    elif opt == 3: # alternate emissions scenario IPCC A1T (low pop. growth, high green energy, low fossil fuel)
        # from https://www.ipcc-data.org/sim/gcm_clim/SRES_TAR/ddc_sres_emissions.html#b1 
        e_GtC_yr = np.array([0, 0,  5.99, 6.9, 8.33, 10, 12.26, 12.6, 12.29, 11.41, 9.91, 8.05, 6.27, 4.31, 0, 0, 0])
    elif opt ==4: # periodic emissions (cos)
        e_GtC_yr = amp*np.cos(((2*np.pi)/period)*(t_yr - 1850)) + e_GtC_yr
    elif opt == 5: # alternate emissions scenario IPCC A1F1 (fossil fuel intensive)
        e_GtC_yr = np.array([0, 0,  5.99, 6.9, 8.65, 11.19, 14.61, 18.66, 23.1, 25.14, 27.12, 29.04, 29.64, 30.32, 0, 0, 0])

    e = np.interp(yr, t_yr, e_GtC_yr)

    return e

def Find_k(F_in, M):
    """Finds values of linear flux constants, k, for fluxes in and out for a 
    given set of fluxes in, out and inital masses. 

    Parameters
    ----------
    F_in: NxM array
        Array with flux-in values
    F_out: NxM array
        Array with flux-out values, positive flux values
    M: M-length vector
        Initial values of mass for each box

    Returns
    -------
    k: NxM array
        Array of k constant values to solve M-box model. k*M gives M ODEs to 
        solve M box model problem
    """
    k_in = F_in/M
    k_out = sum(k_in) # column-wise sum of k values 
    k = k_in - k_out*np.diag(np.ones(len(F_in)))
    return k

def MassFlux(t, M, a, b):
    """Mass flux coupled ODEs for 4-box or 9-box model of carbon cycle in a 
    steady state forced by anthropogenic CO2 emissions

    Parameters
    ----------
    t: t-length array
        Timespan to use for emissions
    M: M-length array
        Initial values of mass for each box
    a: int 
        The forcing scenario to use, int value corresponds to forcing:
            0: No forcings, steady state
            1: IPCC A2 Emissions forcings
            2: IPCC A2 modified so 2000 emissions beyond 2110
            3: some fun periodic forcing
            4: damped sinusoid
    Returns
    -------
    dMdt: M-length array
        Mass flux for one timestep
    """
    global Flux_in_4, Flux_in_9

    ## Flux in and Flux out values, units Gt/yr, to find k values
    if len(M) == 4:
        k = Find_k(Flux_in_4, M0_4)
    else:
        k = Find_k(Flux_in_9, M0_9)

    # The actual matrix multiplication
    dMdt = np.matmul(k, M)

    # Forcings are adding mass to Carbon box in position 0 of the mass vector
    if a==1:
        dMdt[0] += emissions(t) ## adding forcing from emissions function at each timestep
    elif a==2:
        dMdt[0] += emissions(t, opt = 0) ## adding forcing from emissions function at each timestep (maintain 2100 levels)
    elif a==3:
        dMdt[0] += emissions(t, opt = 1) # sine periodic forcing
    elif a==4:
        dMdt[0] += emissions(t, opt = 2) # exponential periodic forcing
    elif a==5:
        dMdt[0] += emissions(t, opt = 3) # try other emissions forcing (A1T)
    elif a==6:
        dMdt[0] += emissions(t, opt = 4) # cos periodic forcing
    elif a==7:
        dMdt[0] += emissions(t, opt = 5) # try other emissions forcing (A1FI)
    

    return dMdt

def rk4(fxy, x0, xf, y0, N):
    """Runge-Kutta integration to solve odes to the 4th order
    Parameters
    ----------
    fxy: function 
        The name of the function containing f(x,y) (e.g. oneode, twoode)
    x0: int OR float 
        Initial values of the independent variable 
    xf: int OR float 
        Final values of the independent variable 
    y0: numpy array 
        Initial value of dependent variable at xo
    N: int 
        Number of intervals to use between xo and xf

    Returns
    -------
    X: numpy array
        Numpy array containing values of the independent variable
    Y: 1D numpy array if one equation solved; M-D numpy array [y1(x) y2(x) ... ] 
            for multiple (M) equations
        The estimated dependent variable at each value of the independent variable
        
    """
    # compute step size and size of output variables
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
    
    # begin computational loop
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
