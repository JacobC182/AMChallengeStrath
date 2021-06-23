#Import libraries/functions
import heyoka as hy
import numpy as np
from numpy import linalg
from numpy.lib.function_base import diff
from numpy.linalg.linalg import norm
from FunctionLibrary import *
from starterODE import ODE
import pykep as pk

GMe = 3.986004407799724e+5

#Choosing AMgrid resolution/spacing
AMspacing = 10


#creating AM-ratio grid
AMgrid = np.linspace(start = 10**-0.5, stop = 10**1.8, num = int((10**1.8 - 10**0.5)/AMspacing), endpoint = True)
AMgrid = np.multiply(AMgrid, 1e-6)
#AMgrid = [10.32583e-3]
#Choose debris file to read from file numbering
fileNum = "001"
#Reading debris state and time vectors
N = int(fileNum)
bigdata = np.loadtxt("data\labels_train.dat")

#Reading debris creation time and state vector
t0 = int(bigdata[N-1][2])*60*60*24
vecInit = pk.par2ic([bigdata[N][3], bigdata[N][4], (bigdata[N][5])*(np.pi/180), (bigdata[N][8])*(np.pi/180), (bigdata[N][7])*(np.pi/180), EccentricAnomalySolver(bigdata[N][6], bigdata[N][4])], GMe)

#Converting creation state vector from 2x3 to 1x6 -> vec0
vec0 = []
for i in range(2):
    for j in range(3):
        vec0.append(vecInit[i][j])
print(vec0)
#Reading debris observation data time vector and state vectors
t, vec = DebrisRead(fileNumber= fileNum)
print(vec)
#Time vector unit conversion days --> secs
t = np.multiply(t,60*60*24)
print(t)

#Creating integrator object
ta = hy.taylor_adaptive(sys = ODE(), state = vec0, high_accuracy = True)#, tol = 1e-10)

#Defining Callable Integrator function that returns the propagated state vectors at the specified time points in the grid
def Integrator(tGrid, AM):
    #resetting integrator initial state and time, +updating AMratio value
    ta.time = t0
    ta.state[:] = vec0
    ta.pars[0] = AM

    #propagating along time grid and returning list of state vectors at grid times
    output = ta.propagate_grid(tGrid)
    print(output)
    return output[4]


#Creating Accuracy list
Accuracy = []

#Looping through Am-ratio grid and propagating orbits for each AM value
for i in range(len(AMgrid)):

    out = Integrator(t, AMgrid[i])
    print(out)
    difference = (np.subtract(abs(out), abs(vec)))

    magdiff = []
    magdiff.append(AMgrid[i]*1e3)

    for j in range(len(difference)):
        
        magdiff.append(np.linalg.norm(difference[j,0:2]))
        magdiff.append(np.linalg.norm(difference[j,3:5]))

    print(magdiff)
    Accuracy.append(magdiff)


print(AMgrid)

np.savetxt(fname = "acctest1.txt", X = Accuracy, delimiter = " ")
