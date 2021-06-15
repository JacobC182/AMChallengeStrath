import time
scriptStartTime = time.time()

#Importing ODE system from ODE function file
from ODE import ODE
#Creating ODE system
Control = [1,1,1,1,1,1,1]
sysODE = ODE(Control)

#Importing Heyoka.py, NumPy, Math Libraries
import heyoka as hy
import numpy as np
import math as ma

#Defining initial conditions - UNITS - KM & KM/s
sysInit = [2.6533e4, 0, 0, 0, 2.2220, 3.1734]
#Creating integrator object
ta = hy.taylor_adaptive(sysODE,sysInit, pars = [0.1], tol = 1e-8)

#printing taylor series object summary
print(ta)

grid = np.linspace(0,946728000,131400) #UNITS - seconds
#propagating along timesteps of grid
#out = ta.propagate_grid(grid)
#print(ta)

#propagating until 30 years
out = ta.propagate_until(t=946728000)
#computing and printing script execution time
scriptExecTime = (time.time() - scriptStartTime)
print('Script Execution time: ' + str(scriptExecTime)[0:6] + ' seconds')