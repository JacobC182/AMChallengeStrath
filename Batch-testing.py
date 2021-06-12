#Andrea Milani Challenge - Main Py Script

#TESTING SCRIPT - Testing the performance impact of batch mode integration

#Importing time module - to measure script execution time
import time
scriptStartTime = time.time()

#Importing ODE system from ODE function file
from ODE import ODE
#Creating ODE system
sysODE = ODE()

#Importing Heyoka.py, NumPy, Math Libraries
import heyoka as hy
import numpy as np
import math as ma

#Defining parameter grid spacing
n = 50

#Creating parameter grid
AMgrid = np.linspace(0.1, 5, n)
AMgrid = np.transpose(AMgrid)

#Defining initial conditions - UNITS - KM & KM/s
#and extrapolating to a constant grid
sysInit = [np.linspace(2.6533e4,2.6533e4,n), np.linspace(0,0,n), np.linspace(0,0,n), np.linspace(0,0,n), np.linspace(2.2220,2.2220,n), np.linspace(3.1734,3.1734,n)]

#Creating integrator object
ta = hy.taylor_adaptive_batch(sysODE, sysInit)

out = ta.propagate_until(np.linspace(200000,200000,n))

#computing and printing script execution time
scriptExecTime = (time.time() - scriptStartTime)
print('Script Execution time: ' + str(scriptExecTime)[0:6] + ' seconds')