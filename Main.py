#---Team Sidereal_Daydreamers - Andrea Milani Challenge 2021---#
#---------------------------------------------------------------
#Importing Libraries/Methods
from mpmath.functions.functions import im
import numpy as np
import heyoka as hy
from numpy.core.function_base import linspace
from FunctionLibrary import CRAMRegressorModel, FileStr, mse, SatRead, orb2rv, rv2orbF
from ODE import ODE
import joblib as jl
from joblib.parallel import delayed
from pykep import ic2par, par2ic

#PART 1 - ESTIMATING THE AREA-TO-MASS RATIO---------------------

Regressor1 = CRAMRegressorModel()     #Creating trained AM-ratio estimator model from function library

initialDebrisState = []     #create empty array to hold the earliest observed state of the debris
initialDebrisTime = []      #create empty array to hold the time corresponding to the debris state
AMratio = []    #create empty AMratio list

for i in range(75,100 +1, 1):
    debrisData = np.loadtxt("data\deb_train\eledebtrain" + FileStr(i) + ".dat")     #Read each debris observation file
    debrisData = np.reshape(debrisData, [-1, 7])

    predictedRatio = Regressor1.predict(debrisData)     #Predicting AM-ratio using Model

    predictedRatio = sum(predictedRatio)/len(predictedRatio)    #Taking average of predictions for multiple observation debris
    AMratio.append(predictedRatio)      #Adding prediction to list

    initialDebrisTime.append(debrisData[0,0])       #adding debris observation time to list
    initialDebrisState.append(debrisData[0,1:7])    #adding debris observed state to list

initialDebrisState = np.reshape(initialDebrisState, [-1,6])     #reshaping initial debris state array from needlessly 3D to 2D

realRatio = np.loadtxt("data\labels_train.dat")[:,1]        #Reading in the REAL ratios from the training data
msError = mse(realRatio[74:100], AMratio)
print("MSE Score: " + str(msError) )

#---------------------------------------------------------------
#PART 2 - PROPAGATING DEBRIS BACK TO SATELLITE------------------

SatData = SatRead()         #Reading all 100 satellite trajectories into 3D array
SatData = np.reshape(SatData, [100,-1,7])
SatTimeGrid = SatData[0,:,0]  #Creating 1D list of satellite trajectory observation time grid
ta = hy.taylor_adaptive(sys = ODE(), state = [0,0,0,0,0,0])     #Creating Heyoka integrator object configured with the dynamical ODE system


#Timestep Callback Function - called at the end of every integration point on the time grid - solves for possible satellite debris causations from satellite trajectories
def CollisionCallback(ta):      #THIS FUNCTION IS GOING TO BE CALLED AT EVERY TIMESTEP FOR EVERY DEBRIS (well over a million times) - Make it good :)
    
    #print("Time: " + str(ta.time/(60*60*24)))
    timeIndex = np.searchsorted(a=SatTimeGrid, v=(ta.time/60*60*24))
    #print(timeIndex)

    debrisState = rv2orbF(ta.state[:])
    satList = SatData[:,timeIndex-1,1:7]
    satList=np.reshape(satList,[100,6])

    for i in range(100):
        if abs(debrisState[0] - satList[i,0]) < 300 and abs(debrisState[1] - satList[i,1]) < 0.1 and abs(debrisState[2] - satList[i,2]) < 1 and abs(debrisState[3] - satList[i,3]) < 1:
            print("Collision Found. " + "Satellite Number " + str(i+1))

    return(True)        #Returns true for no reason other than to not halt the integration process over the time grid - propagate_ callbacks can halt with return(False)


#Numerical Integration Function - A method to propagate over a time grid for each
def Integrator(startTime, startState, AMratio, debrisNum):

    stopTime = -7305 *60*60*24          #Converting from days to seconds, Propagation units KM, KM/s - Treat time variables as (seconds) inside this function
    startTime = startTime *60*60*24

    startVec = orb2rv(startState)       #Converting orbital elements from file to state vector - See Function Library

    ta.time = startTime          #Setting integrator start time (seconds)
    ta.state[:] = startVec       #Setting integratior initial state
    ta.pars[0] = AMratio *1e-6   #Setting ODE parameter (AM-ratio) with units conversion

    while ta.time > stopTime:
        ta.propagate_for(delta_t= -864000, callback = CollisionCallback)   #Back-propagating the debris trajectory for 1 timestep (10 days / 864000s) Until the endTime (day 3645 after EME2000)

    return ta.state[:]


#out = jl.Parallel(n_jobs=-1, prefer="threads")(delayed(Integrator)(initialDebrisTime[i], initialDebrisState[i,:], AMratio[i], i) for i in range(len(AMratio)) )
Integrator(initialDebrisTime[3], initialDebrisState[3,:], AMratio[3], 0)
print(AMratio[3])
#np.savetxt("testOut.txt",out)











#
#TODO
#Create Heyoka non-terminal event function to detect when possible collisions accur
#
#
#Create event function callback method that logs data of possible collisions when they are detected
#
#
#Write implementation of integration in PARALLEL :( (joblib library)
#
#
#
#
#
#Probably lots of other stuff

