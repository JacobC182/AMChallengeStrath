import numpy as np
import heyoka as hy
from FunctionLibrary import *
from ODE import ODE

ta = hy.taylor_adaptive(ODE(), [0,0,0,0,0,0])

propagatedresults = []
propagatedelements = []

observationlist = []

origindata = np.loadtxt("data\labels_train.dat")

for i in range(100):
    print("Debris " + str(i+1) )

    observation = np.loadtxt("data\deb_train\eledebtrain" + FileStr(i+1) + ".dat")
    observation = np.reshape(observation, [-1,7])
    
    ta.time = origindata[i,2]*60*60*24
    ta.state[:] = orb2rv(origindata[i,3:9])
    ta.pars[0] = origindata[i,1]*1e-6

    if len(np.shape(observation)) == 1:     
        ta.propagate_until(t=observation[0]*60*60*24)
        result = ta.state[:]
        propagatedresults.append(result)
    else:
        tGrid = observation[:,0]
        tGrid = np.multiply(tGrid, 60*60*24)
        result = ta.propagate_grid(tGrid)[4]
        for j in range(len(result)):
            propagatedresults.append(result[j])    
        
propagatedresults = np.reshape(propagatedresults, [-1,6])

for i in range(len(propagatedresults)):
    propagatedelements.append(rv2orbF(propagatedresults[i]))

for i in range(100):
    observation = np.loadtxt("data\deb_train\eledebtrain" + FileStr(i+1) + ".dat")
    observation = np.reshape(observation, [-1,7])

    for j in range(len(observation)):
        observationlist.append(observation[j,1:7])

residual = np.subtract(observationlist, propagatedelements)

np.savetxt("residual.txt", abs(residual))