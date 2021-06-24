import time
import heyoka as hy
import numpy as np
import math as ma
import pykep as pk
from joblib.parallel import delayed
from starterODE import ODE
from FunctionLibrary import *
import scipy.optimize as optimization


#Creating ODE system
sysODE = ODE()

#Choosing debris file to read
debrisNumber = "002"

#Reading debris data file
debrisTime, debrisState = DebrisRead(debrisNumber)

#Converting debrisTime from days to seconds
debrisTime = np.multiply(debrisTime, (60*60*24))

#Setting start time as first time debris is spotted
startTime = debrisTime[0]

#setting start state as first state vector from debris observation at start time
startState = debrisState[0][:]

#Removing initial time and state from time and state vectors
debrisTime = np.delete(debrisTime, 0)
for i in range(6):
    debrisState = np.delete(debrisState, [0])


#reshaping state vector from 1D to 2D
debrisState = np.reshape(debrisState, [-1,6])

#extracting positional values only from state vector
debrisStateVel = debrisState[:,3:6]

#creating empty magnitude array
debrisStateMag = []

#calculating magnitude and adding to array, for each given observation
for i in range(len(debrisStateVel)):
    debrisStateMag.append(np.linalg.norm(debrisStateVel[i]))
print(debrisStateMag)

#Creating heyoka integrator object with ODE system
ta = hy.taylor_adaptive(sys = sysODE, state = startState)

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

    #creating empty magnitude storage array
    mag = []

    #calculating magnitude of state vectors for each vector produced by time grid
    for i in range(len(out)):
        mag.append(np.linalg.norm(out[i, 4:6]))

    
    #returning system state vector
    print("Difference:" + str(np.subtract(mag, debrisStateMag)))
    
    return mag

#Setting initial AM-ratio guess
AMguess = 12*1e-6

#Curve fitting ODE function - parameter estimation
optimumRatio, covarianceMatrix = optimization.curve_fit(f = Solution, xdata = debrisTime, ydata = debrisStateMag, p0 = AMguess, bounds = [(10**-0.5)*1e-6, (10**1.8)*1e-6])


print("Fitted AM-Ratio:" + str(optimumRatio[0]*1e6)[0:6])