import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci

# INITIAL VALUES
# initial box masses    
M1 = 725 #atmosphere
M2 = 725 #surface water
M5 = 110 #short-lived biota
M7 = 60 #litter

# fluxes between boxes -> DIFFERENT FROM 'boxmodel.py' TO MAKE NET FLUX 0
F12 = 90
F21 = 90
F15 = 110
F51 = 55
F57 = 55
F71 = 55 #detritus decomposition
F72 = 0

# BOX SYSTEM ODEs
dM1 = [F21, -F12, F51, -F15, F71] #Ffos1
dM2 = [F12, -F21, F72]
dM5 = [F15, -F51, -F57]
dM7 = [F57, -F71, -F72]

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
    return F/M

if __name__ == "__main__":
    M_init = [M1, M2, M5, M7]



