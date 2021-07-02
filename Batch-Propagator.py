import heyoka as hy
import numpy as np
from FunctionLibrary import *
from starterODE import ODE
from joblib.parallel import delayed
from pykep import par2ic

observedDebris = []
propagatedDebris = []
fileNo = 1
ta = hy.taylor_adaptive(sys = ODE(), state = [0,0,0,0,0,0])

for debris in np.loadtxt("data\labels_train.dat")[0:100,:]:
    
    print(str(fileNo) + "% Completed")

    eAnom = EccentricAnomalySolver(debris[6], debris[4])
    r, v = par2ic([debris[3], debris[4], debris[5]*(np.pi/180), debris[8]*(np.pi/180), debris[7]*(np.pi/180), eAnom],3.986004407799724e+5)

    ta.state[0:3] = r
    ta.state[3:7] = v

    ta.time = debris[2] *60*60*24
    ta.pars[0] = debris[1] *1e-6

    fileStr = ""

    if len(str(fileNo)) == 1:
        fileStr = "00" + str(fileNo)
    elif len(str(fileNo)) == 2:
        fileStr = "0" + str(fileNo)
    else:
        fileStr = str(fileNo)

    fileNo += 1

    debrisData = np.loadtxt("data\deb_train\eledebtrain" + fileStr + ".dat")

    observedDebris
    if len(np.shape(debrisData)) == 1:
        debrisData = np.reshape(debrisData, [1,7])
    tGrid = debrisData[:,0]
    tGrid = np.multiply(tGrid, 60*60*24)

    if len(tGrid) == 1:
        ta.propagate_until(t=tGrid)
        propState = ta.state[:]
        propState = np.reshape(propState, [1,6])
    else:
        propState = ta.propagate_grid(tGrid)[4]

    for sighting in debrisData[:,1:7]:
        observedDebris.append(sighting)

    for result in propState:

        orb = rv2orb(result)
        meanAnom = orb[5] - orb[1]*np.sin(orb[5])
        propagatedDebris.append([orb[0], orb[1], orb[2], meanAnom*(180/np.pi), orb[4]*(180/np.pi), orb[3]*(180/np.pi)])


print("saving")
np.savetxt(fname = "observedData.txt", X = observedDebris, delimiter = ",")

np.savetxt(fname = "propagatedData.txt", X = propagatedDebris, delimiter = ",")