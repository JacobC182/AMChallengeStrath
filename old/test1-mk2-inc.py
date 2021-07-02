import time
import heyoka as hy
import numpy as np
from numpy.lib.function_base import cov
from starterODE import ODE
from FunctionLibrary import *
import scipy.optimize as optimization

#Storing script starting time for run-time calculation
ScriptStartTime = time.time()

#Creating ODE system
sysODE = ODE()

#Choosing debris file to read
debrisNumber = "001"

#Reading debris data file into time and state vector lists
debrisTime, debrisState = DebrisRead(debrisNumber)

#Reading debris elements from file into time and orbital elements lists
unusedT, debrisElement = DebrisReadElement(debrisNumber)

#extracting inclination only into list
incList = debrisElement[:,2]

#Converting inclination values from degrees to radians
incList = np.multiply(incList, np.pi/180)

#Removing initial inclination from ecc list
incList = np.delete(incList, 0)

#Converting debrisTime from days to seconds
debrisTime = np.multiply(debrisTime, (60*60*24))

#Setting start time as first time debris is spotted
startTime = debrisTime[0]

#setting start state as first state vector from debris observation at start time
startState = debrisState[0][:]

#Removing initial time and state from time and state vectors
debrisTime = np.delete(debrisTime, 0)


#Creating heyoka integrator object with ODE system
ta = hy.taylor_adaptive(sys = sysODE, state = startState, tol=1e-20)

#defining ODE result returning function
def Solution(t, AM):
    #resetting integrator state and time
    ta.state[:] = startState
    ta.time = startTime
    ta.pars[0] = AM

    #propagating over time points
    out = ta.propagate_grid(t)

    #removing unnecessary data from integrator output
    out = out[4]

    #creating empty inclination storage array
    inc = []

    #calculating inclination for each time-grid state vector
    for i in range(len(out)):
        inc.append(rv2orb(out[i])[2])
    
    #printing ON-THE-FLY results
    print("Am-Ratio:  " + str(AM*1e6)[0:6] )
    print("Inclination Error: " + str(np.subtract(incList, inc)) )
    print("<---------------------------->")
    
    #returning inclination (ydata)
    return inc

#Setting initial AM-ratio guess
AMguess = 4*1e-6

#Curve fitting ODE function - parameter estimation
optimumRatio, covarianceMatrix = optimization.curve_fit(f = Solution, xdata = debrisTime, ydata = incList, p0 = AMguess, bounds = [(10**-0.5)*1e-6, (10**1.8)*1e-6])
#optimumRatio, covarianceMatrix = optimization.curve_fit(f = Solution, xdata = debrisTime, ydata = eccList, bounds = [(10**-0.5)*1e-6, (10**1.8)*1e-6])
#Printing final result
print("Fitted AM-Ratio:" + str(optimumRatio[0]*1e6))

#Printing script run-time
print("Script Finished In: " + str(time.time() - ScriptStartTime)[0:7] + "s")
