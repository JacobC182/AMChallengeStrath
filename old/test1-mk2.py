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
debrisNumber = "001"

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
print(debrisState)

#reshaping state vector from 1D to 2D
debrisState = np.reshape(debrisState, [-1,6])


ta = hy.taylor_adaptive(sys = sysODE, state = startState)

#defining ODE result returning function
def Solution(t, AM):
    #resetting integrator state and time
    ta.state[:] = startState
    ta.time = startTime
    ta.pars[0] = AM
    #propagating over time points
    out = ta.propagate_grid(t)
    #returning system state vector
    return out[4]

#Setting initial AM-ratio guess
AMguess = 10*1e-6

#Curve fitting ODE function - parameter estimation
optimumRatio, covarianceMatrix = optimization.curve_fit(f = Solution, xdata = debrisTime, ydata = debrisState, p0 = AMguess, bounds = [(10**-0.5)*1e-6, (10**1.8)*1e-6])

print(optimumRatio, covarianceMatrix)