#Andrea Milani Challenge - Main Py Script

#Importing time module - to measure script execution time
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
ta = hy.taylor_adaptive(sysODE,sysInit, pars = [6], tol = 1e-7)

#printing taylor series object summary
print(ta)
print("Taylor Decomposition Size: " + str(len(ta.decomposition)-12) )

#Calculating Orbital Period
orbRadius = (sysInit[0]**2 + sysInit[1]**2 + sysInit[2]**2)**0.5
GME = np.double(3.986004407799724e+5)   #GM-Earth (KM3/s2)
orbPeriod = 2*ma.pi*ma.sqrt((orbRadius**3)/GME)

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


#Printing Initial orbital period - testing
print("Initial Orbital Period: " + str(orbPeriod*0.00001157407)[0:6] + " Days")
print("Propagation Period: " + str(grid[-1]*0.00001157407)[0:6] + " Days")
print("No. Of Orbits:" + str(grid[-1]/orbPeriod)[0:7])


#saving propagator output to txt file
np.savetxt("propOut.txt", out[4],delimiter=',')

#saving time grid to txt file
np.savetxt("timeOut.txt", grid[:])

#calculating magnitude of body position
posMagnitude = np.sqrt((out[4][:, 3]**2) + (out[4][:, 4]**2) + (out[4][:, 5]**2))

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