import numpy as np

## Defining functions, 
def emissions(yr):
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
    # changed last 3 variables in vector to the 2000 year emission value (previously all 0)
    e_GtC_yr = np.array([0, 0,  6.875, 8.125, 9.375, 12.5, 14.375, 16.25, 17.5, 19.75, 21.25, 23.125, 26.25, 28.75, 8.125, 8.125, 8.125])

    e = np.interp(yr, t_yr, e_GtC_yr)

    return e

def MassFlux_forced(t, M):
    """
    Mass flux coupled ODEs for 4-box model of carbon cycle in a steady state
    """
    # add emissions to box 1 
    M[0] += emissions(t)
    ## Flux in and Flux out values, units Gt/yr
    Flux_in = np.array([[0, F21, F51, F71], \
                    [F12, 0, F52, F72], \
                    [F15, F25, 0, F75],\
                    [F17, F27, F57, 0]])
    Flux_out = np.array([[0, F12, F15, F17], \
                    [F21, 0, F25, F27], \
                    [F51, F52, 0, F57],\
                    [F71, F72, F75, 0]])
    Mass_Flux_total = sum(Flux_in - Flux_out)

    ## find K values for intial Steady State fluxes given
    k_in = Flux_in.T/M
    k_out = sum(Flux_out/M)
    k = k_in.T - k_out*np.diag(np.ones(4))

    ## Matrix multiplication to create linear system of equations
    dMdt = np.matmul(k, M) # this result is the steady state solution
    return dMdt

def MassFlux_Steady(t, M):
    """
    Mass flux coupled ODEs for 4-box model of carbon cycle in a steady state
    """
    ## Flux in and Flux out values, units Gt/yr
    Flux_in = np.array([[0, F21, F51, F71], \
                    [F12, 0, F52, F72], \
                    [F15, F25, 0, F75],\
                    [F17, F27, F57, 0]])
    Flux_out = np.array([[0, F12, F15, F17], \
                    [F21, 0, F25, F27], \
                    [F51, F52, 0, F57],\
                    [F71, F72, F75, 0]])
    Mass_Flux_total = sum(Flux_in - Flux_out)

    ## find K values for intial Steady State fluxes given
    k_in = Flux_in.T/M
    k_out = sum(Flux_out/M)
    k = k_in.T - k_out*np.diag(np.ones(4))

    ## Matrix multiplication to create linear system of equations
    dMdt = np.matmul(k, M) # this result is the steady state solution
    return dMdt


def FindK(F ,t):
    """
    Function to integrate given fluxes with time to see if system is in steady state. It is already, makes sense no feed back to M in equations 
    (can't be done without k anyway))
    """
    dMdt = F 
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
