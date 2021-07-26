#---Team Sidereal_Daydreamers - Andrea Milani Challenge 2021---#
#---------------------------------------------------------------
#Importing Libraries/Methods
import numpy as np
import heyoka as hy
from numpy.core.function_base import linspace
from FunctionLibrary import CRAMRegressorModel, FileStr, mse, SatRead, orb2rv, rv2orbF, XrayVision
from ODE import ODE
import joblib as jl
from joblib.parallel import delayed
from pykep import ic2par, par2ic
import time
#PART 1 - OUTLIER DETECTION // REMOVING NOISY OBSERVATIONS------

DebrisDetector = XrayVision()       #Creating trained Outlier detection model from function library
debData = []

for i in range(1,25 +1, 1):
    debrisData = np.loadtxt("data\deb_train\eledebtrain" + FileStr(i) + ".dat")     #Read each debris observation file

    if len(np.shape(debrisData)) == 1:      #Appending debris data and AM ratio to input and output training arrays respectivel
            debData.append(debrisData)
    else:
        for j in range(len(debrisData)):        #Multiple observations from the same debris are given the same true AM-ratio from that debris in the 2D training array format
            debData.append(debrisData[j,:])

noiseResults = DebrisDetector.predict(debData)
print(noiseResults)
#---------------------------------------------------------------
#PART 2 - ESTIMATING THE AREA-TO-MASS RATIO---------------------

#TIMING
st1 = time.time()

Regressor1 = CRAMRegressorModel()     #Creating trained AM-ratio estimator model from function library

initialDebrisState = []     #create empty array to hold the earliest observed state of the debris
initialDebrisTime = []      #create empty array to hold the time corresponding to the debris state
AMratio = []    #create empty AMratio list

for i in range(1,25 +1, 1):
    debrisData = np.loadtxt("data\deb_train\eledebtrain" + FileStr(i) + ".dat")     #Read each debris observation file
    debrisData = np.reshape(debrisData, [-1, 7])

    predictedRatio = Regressor1.predict(debrisData)     #Predicting AM-ratio using Model

    predictedRatio = sum(predictedRatio)/len(predictedRatio)    #Taking average of predictions for multiple observation debris
    AMratio.append(predictedRatio)      #Adding prediction to list

    initialDebrisTime.append(debrisData[0,0])       #adding debris observation time to list
    initialDebrisState.append(debrisData[0,1:7])    #adding debris observed state to list

initialDebrisState = np.reshape(initialDebrisState, [-1,6])     #reshaping initial debris state array from needlessly 3D to 2D

realRatio = np.loadtxt("data\labels_train.dat")[:,1]        #Reading in the REAL ratios from the training data
msError = mse(realRatio[0:25], AMratio)
print("MSE Score: " + str(msError) )

#---------------------------------------------------------------
#PART 3 - PROPAGATING DEBRIS BACK TO SATELLITE------------------

SatData = SatRead()         #Reading all 100 satellite trajectories into 3D array
SatData = np.reshape(SatData, [100,-1,7])


SatTimeGrid = SatData[0,:,0]                                #Creating 1D list of satellite trajectory observation time grid
SatTimeGrid = np.multiply(SatTimeGrid, 60*60*24)            #Converting time grid from days to seconds


ta = hy.taylor_adaptive(sys = ODE(), state = [0,0,0,0,0,0])     #Creating Heyoka integrator object configured with the dynamical ODE system



#Numerical Integration Function - A method to propagate over a time grid for each debris
def Integrator(startTime, startState, AMratio, debrisNum):         
    
    startTime = startTime *60*60*24     #Converting from days to seconds, Propagation units KM, KM/s - Treat time variables as (seconds) inside this function

    timeIndex = np.searchsorted(a=SatTimeGrid, v=startTime)     #Getting array index number for start time in time list
    tGrid = SatTimeGrid[0:timeIndex]                            #Creating time grid from start time back to earliest trajectory time
    tGrid = np.flip(tGrid)                                      #Flipping time grid to go back-in-time (back-propagation)

    sat = SatData[:,0:timeIndex,1:7]
    sat = np.flip(sat, axis=1)
    
    startVec = orb2rv(startState)       #Converting orbital elements from file to state vector - See Function Library
    #print(DeNoise.predict([startState]))
    #startVec =  orb2rv( DeNoise.predict([startState])[0] )         #DENOISING INPUT OBSERVATION

    ta.time = startTime          #Setting integrator start time (seconds)
    ta.state[:] = startVec       #Setting integratior initial state
    ta.pars[0] = AMratio *1e-6   #Setting ODE parameter (AM-ratio) with units conversion

    vecOut = ta.propagate_grid(grid = np.flip(SatTimeGrid))[4]  #Back-propagating the debris trajectory over the time grid (timestep = 10 days / 864000s) Until the endTime (day 3645 after EME2000)

    out = []

    for vec in vecOut:
        out.append(rv2orbF(vec))

    out = np.reshape(out, [-1,6])

    nCollisions = 0

    #SETTING TOLERANCE VALUES
    kmTol = 700
    eccTol = 0.2
    incTol = 10
    degTol = 20

    for i in range(len(tGrid)):         #Everything inside this loop has to be VERY FAST - Make it good :)
        for j in range(100):
            if abs(out[i,0] - sat[j,i,0]) < kmTol and abs(out[i,1] - sat[j,i,1]) < eccTol and abs(out[i,2] - sat[j,i,2]) < incTol and abs(out[i,3] - sat[j,i,3]) < degTol and abs(out[i,4] - sat[j,i,4]) < degTol and abs(out[i,5] - sat[j,i,5]) < degTol:
                print("Collision Found: Sat-" + str(j+1) + "  Time: " + str(tGrid[i] /(60*60*24)) + " days")
                nCollisions += 1

    return(print("No. of Detections: " + str(nCollisions) + "  For debris: " + str(debrisNum +75)))


#out = jl.Parallel(n_jobs=-1, prefer="threads")(delayed(Integrator)(initialDebrisTime[i], initialDebrisState[i,:], AMratio[i], i) for i in range(len(AMratio)) )
debNumber = 75  -75#Debris chosen number

#Integrator(initialDebrisTime[debNumber], initialDebrisState[debNumber,:], AMratio[debNumber], debNumber)


print("AM-Ratio: " + str(AMratio[debNumber])[0:6])
#np.savetxt("testOut.txt",out)


print("Run Time: " + str(time.time()- st1)[0:7] + "s")








#---------------------------------------------------------------
#TODO
#Create collision detection system
#DONE
#
#
#
#
#Write implementation of integration in PARALLEL :( (joblib library)
#DONE
#
#
#
#
#DE-NOISE debris observation data :()
#
#
#
#
#Probably lots of other stuff

