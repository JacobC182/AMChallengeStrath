#Andrea Milani Challenge - Main Py Script

#Import Libraries
import time
import heyoka as hy
import numpy as np
import math as ma
import pykep as pk
from joblib.parallel import delayed
from starterODE import ODE


scriptStartTime = time.time()
#Creating ODE system
sysODE = ODE()

#

for i in range(100):

    fileNo = i
    if len(str(fileNo)) == 1:
        fileNo = int("00" + str(fileNo))
    elif len(str(fileNo)) == 2:
        fileNo = int("0" + str(fileNo))

    debFile = open("\data\deb_train\eledebtrain" + str(fileNo) + ".dat", "r")

    debFileLength = len(debFile.readlines())

    

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