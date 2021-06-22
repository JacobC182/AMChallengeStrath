#Import libraries/functions
import heyoka as hy
import numpy as np
from numpy.lib.function_base import diff
from FunctionLibrary import *
from starterODE import ODE

#Choosing AMgrid resolution/spacing
AMspacing = 10


#creating AM-ratio grid
AMgrid = np.linspace(start = 10**-0.5, stop = 10**1.8, num = int((10**1.8 - 10**0.5)/AMspacing), endpoint = True)


#Choose debris file to read from file numbering
fileNum = "001"
#Reading debris state and time vectors
t, vec = DebrisRead(fileNumber= fileNum)

t2 = t
t2.pop(0)
vec2 = np.delete(vec, 0,0)


#Creating integrator object
ta = hy.taylor_adaptive(sys = ODE(), state = vec[0], tol = 1e-10)

#Defining Callable Integrator function that returns the propagated state vectors at the specified time points in the grid
def Integrator(tGrid, AM):

    ta.time = tGrid[0]
    ta.state[:] = vec[0]
    ta.pars[0] = AM

    return ta.propagate_grid(tGrid)[4]


#Creating Accuracy list
Accuracy = []

#Looping through Am-ratio grid and propagating orbits for each AM value
for i in range(len(AMgrid)):

    out = Integrator(t2, AMgrid[i])

    rd = (np.subtract(abs(out), abs(vec2)))

    print(rd)

    #Accuracy.append(difference)
print(AMgrid)
#Accuracy = np.reshape(Accuracy, [len(t2), 6])
#print(Accuracy)
#np.savetxt(fname = "acctest1.txt", X = Accuracy, delimiter = " ")