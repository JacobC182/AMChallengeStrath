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
    debVector.append(pk.par2ic([debElements[i][1], debElements[i][2], debElements[i][3], debElements[i][6], debElements[i][5], eAnom]))


def Solution(t, initial, AM):

    ta = hy.taylor_adaptive(sys = sysODE, state = initial, pars = [AM], tol=1e-8)

    stepSize = 600

    nPoints = int(abs(timeVector[-1] - timeVector[0]) / stepSize)
    grid = np.linspace(start = timeVector[0], stop = timeVector[-1], num = nPoints, endpoint=True)

    grid = timeVector[1:-1]

    output = ta.propagate_grid(grid)

    return output[4]

def realSolution(t, AM):

    initial = []

    for i in range(2):
        for j in range(3):
            initial.append(debVector[0][i][j])

    print(initial)
    return(Solution(t, list(initial), AM))

AMguess = 1


initial0 = []
for k in range(3):
    for i in range(2):
        for j in range(3):
            initial0.append(debVector[k][i][j])

initial0 = np.reshape(initial0, [3, 6])

optimumRatio, covarianceMatrix = optimization.curve_fit(f = realSolution, xdata = timeVector[1:-1], ydata = initial0, p0 = AMguess, bounds = [10**-0.5, 10**1.8])

print(optimumRatio)