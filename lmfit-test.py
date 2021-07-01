#TESTING - doesnt yield any good results

import time
import heyoka as hy
import numpy as np
from numpy.lib.function_base import cov
from starterODE import ODE
from FunctionLibrary import *
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.optimize import least_squares

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

#extracting eccentricity only into list
eccList = debrisElement[:,1]

#Removing initial eccentricity from ecc list
eccList = np.delete(eccList, 0)

#extracting inclination only into list
incList = debrisElement[:,2]

#Removing initial inclination from inc list
incList = np.delete(incList, 0)

#Converting debrisTime from days to seconds
debrisTime = np.multiply(debrisTime, (60*60*24))

#Setting start time as first time debris is spotted
startTime = debrisTime[0]

#setting start state as first state vector from debris observation at start time
startState = debrisState[0][:]

#Removing initial time and state from time and state vectors
debrisTime = np.delete(debrisTime, 0)
debrisState = np.delete(debrisState, 0, axis=0)

#Creating heyoka integrator object with ODE system
ta = hy.taylor_adaptive(sys = sysODE, state = startState, tol=1e-8)

tGrid = debrisTime


#defining ODE result returning function
def Solution(AM):
    #resetting integrator state and time
    ta.state[:] = startState
    ta.time = startTime
    ta.pars[0] = AM

    #propagating over time points
    out = ta.propagate_grid(tGrid)

    #removing unnecessary data from integrator output
    out = out[4]

    #creating empty eccentricity storage array
    ecc = []
    inc = []
    #calculating eccentricity for each time-grid state vector
    for i in range(len(out)):
        ecc.append(rv2orb(out[i])[1])
        inc.append(rv2orb(out[i])[2])

    residual = [np.subtract(eccList, ecc), np.subtract(incList, inc)]
    residual = np.subtract(out, debrisState)
    #printing ON-THE-FLY results
    print("Am-Ratio:  " + str(AM*1e6)[0:6] )
    print("Error: " + str(residual) )
    print("<---------------------------->")

    #Converting residual output from Python list to NP array
    residualArray = np.array(residual)
 
    #returning eccentricity (ydata)
    return residualArray.flatten()

#Setting initial AM-ratio guess
AMguess = 10*1e-6


#Minimizer-residual fitting
Result = least_squares(fun = Solution, x0 = AMguess, bounds=((10**-0.5)*1e-6, (10**1.8)*1e-6) )



#Printing final result
#print("Fitted AM-Ratio:" + str(optimumRatio[0]*1e6))

#Printing script run-time
print("Script Finished In: " + str(time.time() - ScriptStartTime)[0:7] + "s")
