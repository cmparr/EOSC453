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
FOSSILFUELS = 10000 # 10?

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
F78 = 3
F79 = 1
F81 = 3
F91 = 1
# Ffos1 = ?

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

def boxODE(b):
    # how to turn constructing the ode into one function
    return b

if __name__ == "__main__":
    M = [M1, M2, M3, M4, M5, M6, M7, M8, M9]
