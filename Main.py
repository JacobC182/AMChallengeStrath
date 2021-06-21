#Andrea Milani Challenge - Main Py Script

#Importing time module - to measure script execution time
import time
scriptStartTime = time.time()

#Importing ODE system from ODE function file
from starterODE import ODE
#Creating ODE system

sysODE = ODE()

#Importing Heyoka.py, NumPy, Math Libraries
import heyoka as hy
import numpy as np
import math as ma

#Defining initial conditions - UNITS - KM & KM/s
sysInit = [2.6533e4, 0, 0, 0, 2.2220, 3.1734]


#Creating integrator object
ta = hy.taylor_adaptive(sysODE,sysInit, pars = [6], tol = 1e-7)

#printing taylor series object summary
print(ta)
print("Taylor Decomposition Size: " + str(len(ta.decomposition)-12) )

#Calculating Orbital Period

#creating timestep grid for integrator - UNITS - seconds
#user controlled start/end/step variables below
startTime = 0 #EME2000 Datum
endTime = 946728000
stepSize = 1

nPoints = int((endTime - startTime) / stepSize)
#creating linspace timestep grid
grid = np.linspace(start=startTime, stop=endTime, num=nPoints)

#propagating along timesteps of grid
out = ta.propagate_grid(grid)
#print(ta)

#computing and printing script execution time
scriptExecTime = (time.time() - scriptStartTime)
print('Script Execution time: ' + str(scriptExecTime)[0:6] + ' seconds')


#importing motplotlib library pyplot method (MATLAB-Like Plotting Module)
import matplotlib.pyplot as plt

#3d positional orbit plot
fig = plt.figure(figsize=(9, 9))
ax = plt.axes(projection ='3d')

#plotting body position
plt.plot(out[4][:, 3], out[4][:, 4],out[4][:, 5])

#enable gridlines
plt.grid()
#show plot figure window
plt.show()