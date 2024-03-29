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

#Reading debris data file
debElements = np.loadtxt("data\deb_train\eledebtrain001.dat")

#create empty debris state vector and time arrays
debVector = []
timeVector = []

for i in range(len(debElements)):

    eAnom = EccentricAnomalySolver(debElements[i][4], debElements[i][2])
    timeVector.append(debElements[i][0])
    debVector.append(pk.par2ic([debElements[i][1], debElements[i][2], debElements[i][3]*(np.pi/180), debElements[i][6]*(np.pi/180), debElements[i][5]*(np.pi/180), eAnom]))

timeVector = np.multiply(timeVector, (24*60*60))
timeVector0 = np.delete(timeVector, 0)
print(timeVector0)

def Solution(t, initial, AM):

    ta = hy.taylor_adaptive(sys = sysODE, state = initial, pars = [AM], tol=1e-8)

    stepSize = 600

    nPoints = int(abs(timeVector[-1] - timeVector[0]) / stepSize)
    grid = np.linspace(start = timeVector[0], stop = timeVector[-1], num = nPoints, endpoint=True)

    grid = timeVector0
  
    output = ta.propagate_grid(grid)
    #output = ta.propagate_until(t)
    print(output[4])
    print("-------------------")
    return output[4]

def realSolution(t, AM):

    initial = []

    for i in range(2):
        for j in range(3):
            initial.append(debVector[0][i][j])

    #print(initial)
    return(Solution(t, initial, AM))

AMguess = 10*1e-6


initial0 = []
for k in range(3):
    for i in range(2):
        for j in range(3):
            initial0.append(debVector[k+1][i][j])

initial0 = np.reshape(initial0, [3, 6])

initial0 = initial0[:,1]
print(initial0)



optimumRatio, covarianceMatrix = optimization.curve_fit(f = realSolution, xdata = timeVector0, ydata = initial0[:][1], p0 = AMguess, bounds = [(10**-0.5)*1e-6, (10**1.8)*1e-6])

print(optimumRatio, covarianceMatrix)