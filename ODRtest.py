import numpy as np
import heyoka as hy
from lmfit import Model
import matplotlib.pyplot as plt
from FunctionLibrary import *
import pykep as pk
#Importing ODE creator function from file
from starterODE import ODE

#creating ODE system
sysODE = ODE()

#Reading debris data file
debElements = np.loadtxt("data\deb_train\eledebtrain001.dat")

#create empty debris state vector and time arrays
debVector = []
timeVector = []

#Creating initial A/M ratio guess variable
AM0 = 10

#Parsing/creating vectors of time and state vector arrays
for i in range(len(debElements)):

    eAnom = EccentricAnomalySolver(debElements[i][4], debElements[i][2])
    timeVector.append(debElements[i][0])
    debVector.append(pk.par2ic([debElements[i][1], debElements[i][2], debElements[i][3], debElements[i][6], debElements[i][5], eAnom]))

#Creating initial state vector variable
initVector = np.reshape(debVector[0], 6)

debVector = np.reshape(debVector, [4,6])

#Creating Heyoka integrator object
ta = hy.taylor_adaptive(sys = sysODE, state = initVector, pars = [AM0], tol = 1e-9)

#Defining ODE solution function
def Solution(AM):
    
    ta.time = timeVector[0]
    ta.state = initVector
    ta.pars[0] = [AM]

    return ta.propagate_grid(timeVector)[4]


ODEmodel = Model(Solution, param_names="AM")

ODEmodel.make_params(AM=10)

result = ODEmodel.fit(debVector)
