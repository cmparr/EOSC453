import numpy as np
import matplotlib.pyplot as plt

# CONSTANTS
#time span in years
yr_MIN = 1800
yr_MAX = 2150


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


if __name__ == "__main__":
    #the function linspace makes a linearly spaced vector between the min and max years of n points (n=200 below)
    yr = np.linspace(yr_MIN, yr_MAX, 200)

    #call the function emissions
    e = emissions(yr)

    plt.plot(yr, e)
    plt.xlabel('Time (years')
    plt.ylabel('CO_2 emissions (GtC/yr)')
    plt.title('Modified IPPC scenario A2')

    plt.savefig('IPPC_A2_emissions.pdf')
    plt.close()
