import heyoka as hy
import numpy as np
from FunctionLibrary import *
from starterODE import ODE

#Choose debris number from data ------------------------!
orbitNumber = "001"

#Reading in time and state vector lists from file
timeVector, stateVector = DebrisRead(orbitNumber)

#Reading time and orbital elements list from file
unusedT, elementList = DebrisReadElement(orbitNumber)

#Reading true/correct AM-ratio from file
trueAM = DebrisLabel(orbitNumber)[1]

#Extracting starting time
startTime = timeVector[0]
startTime = startTime*60*60*24

#Extracting starting state vector
startState = stateVector[0]

#Extracting ending time
endTime = timeVector[1]
endTime = endTime*60*60*24

#Extracting correct eccentricity value
trueInc = elementList[1,2]
trueInc = trueInc * np.pi/180

#Creating ODE system object from ODE file
ODE = ODE()

#Creating Heyoka integrator object
ta = hy.taylor_adaptive(sys = ODE, state = [0,0,0,0,0,0])

#Defining callable integration function that propagates for any AM ratio given
def Solve(AM):

    ta.state[:] = startState
    ta.time = startTime
    ta.pars[0] = AM

    ta.propagate_until(t=endTime)

    ElementResult = rv2orb(ta.state)

    return ElementResult[2]

AMlist = np.linspace((10**-0.5)*1e-6, (10**1.8)*1e-6, 200, endpoint=True)

results = []

for i in AMlist:
    results.append(abs(Solve(i)-trueInc))

import matplotlib.pyplot as plt

plt.plot([trueAM, trueAM], [0, 1])
plt.plot(np.multiply(AMlist, 1e6), results)

print("Global Minimum (Error):")
print(min(results))
print("With AM-Ratio:")
print(AMlist[results.index(min(results))]*1e6)
print("True AM-Ratio:")
print(trueAM)

plt.show()