import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci

# INITIAL VALUES
# initial box masses    
M1 = 725 #atmosphere
M2 = 725 #surface water
M3 = 3 #surface biota
M4 = 37675 #intermediate and deep water
M5 = 110 #short-lived biota
M6 = 450 #long-lived biota
M7 = 60 #litter
M8 = 1350 #soil
M9 = 160 #peat
MFOSSILFUELS = 10000 # 10?

# fluxes between boxes
F12 = 89
F21 = 90
F23 = 40
F32 = 36
F34 = 4 #detritus
F42 = 42
F24 = 38
F15 = 110
F51 = 55
F61 = 0 #deforestation
F56 = 15
F57 = 40
F67 = 15
F71 = 50 #detritus decomposition
F72 = 1
F78 = 3
F79 = 1
F81 = 3
F91 = 1
# Ffos1 = ?

# BOX SYSTEM ODEs
dM1 = [F21, -F12, F51, -F15, F61, F71, F81, F91] #+ Ffos1
dM2 = [F12, -F21, F42, -F24, F32, -F23, F72]
dM3 = [F23, -F32, -F34]
dM4 = [F34, -F42, F24]
dM5 = [F15, -F51, -F56, -F57]
dM6 = [F56, -F67, -F61]
dM7 = [F57, -F71, -F72, -F78, -F79, F67]
dM8 = [F78, -F81]
dM9 = [F79, -F91]

# FUNCTION DEFINITIONS
def linflux(k, M):
    """Calculates the box flux (assuming a linear flux law)

    Parameters
    ----------
    k: int
        rate coefficient
    M: int
        mass of box (that mass is moving from)

    Returns
    -------
    k*m: int
        flux between boxes (F)

    """
    return k*M

def ratecoeff(F,M):
    """Calculates the box flux (assuming a linear flux law)

    Parameters
    ----------
    F: int
        flux between boxes
    M: int
        mass of box (that mass is moving from)

    Returns
    -------
    F/M: int
        rate coefficient (k)

    """
    return np.divide(F,M)

# def boxODE(b):
#     # how to turn constructing the ode into one function
#     return b

if __name__ == "__main__":
    # intial mass vector
    M = [M1, M2, M3, M4, M5, M6, M7, M8, M9]

    # initial time
    t0 = 0

    # max integration time
    tmax = 1000

    # time step
    dt = 1 # yearly?

    # initial y
    dMdt = dM1 # TODO put this in a loop to go over all dM/dt eq.
    y0 = np.sum(dMdt)

    # mass of current box 
    M_step = M1 # TODO initiate mass in loop M[n]

    #initialize
    steps = int(np.round(tmax/dt))
    y = np.empty((steps+1,))
    t = np.empty((steps+1,))
    y[0] = y0 # starts at 0
    t[0] = t0
    
    #march forward in time
    for step in range(steps):
        # initialize vector to store info for each flux in dM/dt
        ks = np.zeros(len(dMdt))
        Ms = np.zeros(len(dMdt))

        for n in range(len(dMdt)):
            F_step = dMdt[n] # current flux
            k_step = ratecoeff(F_step, M_step) # calculate rate coefficient using previous F and M values

            ks[n] = k_step
            Ms[n] = M_step

        dMdt = sci.odeint(linflux, y[-1], t, args=(ks, Ms))
        y[step+1] = np.sum(dMdt)
        t[step+1] = t[step] + dt  

    plt.plot(t,y)
    plt.xlabel('time')
    plt.ylabel('dM/dt')

    plt.show()
        



