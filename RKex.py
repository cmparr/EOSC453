'''   
 This script and the accompanying functions is designed to help you learn about the numerical solution of initial value problems (ODEs).

The algorithm is as follows:

1. Define Equation(s) to be solved. This will be a 1st order ODE or system of ODEs. 
2. Define or input appropriate control parameters (such as the length of time for the integration and the tolerance for convergence).
3. Define the vector of initial conditions -- value(s) of the function at time = 0
4. Call the ODE solver
5. Plot your results

Note that typing '<function>??' at the command line is very usefl!

A global variable can be passed to all called functions; 'a' and 'b' are constants in the ODE in the function 'oneode.ipynb': dy/dt = -a*y^b  

Note: if b = 1 then this ode is linear and the analytical solution can be obtained by separation of variables (this is raddioactive decay): y(t) =exp(-at) 

The larger b is the more nonlinear the equations to be solved are. In the coupled ode example the nonlinear coupling has a large effect-- 
Increase or decrease b from 1 with small steps (< than 1) and explore the behavior; You may have to increase the number of time steps 
below to achieve resolution.
'''

import numpy as np
import matplotlib.pyplot as plt

global a, b

a = -1 #constant
b = 0.8

def rk4(fxy, x0, xf, y0, N):
    
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
    M = np.max(np.shape(y0))
    Y = np.zeros((N+1, M))*1j #make complex by multiplying by 1j; this way can add complex values to this during integration

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

def oneode(t,y):
    
    # input:
    #     t: float
    #     y: numpy array with one value

    # output:
    #     ydot: numpy array of same size as y
    
    global a, b #we are bringing a defined in the main program.
    ydot = np.array([a * y[0]**b])
    
    return ydot

def twoodes(t,y):
    
    # input:
    #     t: float
    #     y: numpy array with length = number of equations being solved

    # output:
    #     ydot: numpy array of same size as y
    #            made of complex values because np.power(negative, fraction) returns nan, but np.power(negative, fraction, dtype = np.complex128) returns a complex value

    global a, b #we are bringing a defined in the main program.
    ydot = np.array([a*np.power(y[1],b, dtype = np.complex128), -a*np.power(y[0],b, dtype = np.complex128)]);
    
    return ydot

#time interval for integration
time_min = 0
time_max = 10

#define a linearly-spaced vector with n points
n = 1000
timespan = -np.linspace(time_min, time_max, n)

#array of initial conditions defined in the function 'oneode.ipynb'
y0 = np.array([100.])

#array of initial conditions defined in the function 'twoodes.ipynb'
yy0 = np.array([100.,100.])


[t,y] = rk4(fxy = oneode,
            x0 =  time_min,
            xf = time_max,
            y0 = y0,
            N = n)

[tt,yy] = rk4(fxy = twoodes,
            x0 =  -10,
            xf = 0,
            y0 = yy0,
            N = n)


plt.plot(t,y)
plt.xlabel('Time')
plt.ylabel('y')
plt.title('One ODE Solution')

plt.show()

plt.plot(tt,yy[:,0], label = 'y1')
plt.plot(tt,yy[:,1], label = 'y2')
plt.plot(tt, yy[:,0] + yy[:,1], label = 'sum')
plt.xlabel('Time')
plt.ylabel('y')
plt.title('Coupled ODE Solution')
plt.legend()

plt.show()
