#---Team Sidereal_Daydreamers - Andrea Milani Challenge 2021---#
#---------------------------------------------------------------
#Importing Libraries/Methods
from mpmath.functions.functions import im
import numpy as np
import heyoka as hy
from FunctionLibrary import *
from ODE import ODE

#PART 1 - ESTIMATING THE AREA-TO-MASS RATIO---------------------

Regressor1 = CRAMForestRegressorModel()     #Creating trained AM-ratio estimator model from function library

initialDebrisState = []     #create empty array to hold the earliest observed state of the debris
initialDebrisTime = []      #create empty array to hold the time corresponding to the debris state
AMratio = []    #create empty AMratio list

for i in range(1,100 +1, 1):
    debrisData = np.loadtxt("data\deb_train\eledebtrain" + FileStr(i) + ".dat")     #Read each debris observation file
    debrisData = np.reshape(debrisData, [-1, 7])

    predictedRatio = Regressor1.predict(debrisData)     #Predicting AM-ratio using Model

    predictedRatio = sum(predictedRatio)/len(predictedRatio)    #Taking average of predictions for multiple observation debris
    AMratio.append(predictedRatio)      #Adding prediction to list

    initialDebrisTime.append(debrisData[0,0])       #adding debris observation time to list
    initialDebrisState.append(debrisData[0,1:7])    #adding debris observed state to list

initialDebrisState = np.reshape(initialDebrisState, [-1,6])     #reshaping initial debris state array from needlessly 3D to 2D
#---------------------------------------------------------------
#PART 2 - PROPAGATING DEBRIS BACK TO SATELLITE------------------

ta = hy.taylor_adaptive(sys = ODE(), state = [0,0,0,0,0,0])     #Creating Heyoka integrator object configured with the dynamical ODE system

def Integrator(startTime, startState, AMratio):

    ta.time = startTime *60*60*24   #Setting integrator start time (seconds)
    ta.state[:] = startState        #Setting integratior initial state
    ta.pars[0] = AMratio *1e-6      #Setting ODE parameter (AM-ratio) with units conversion

print(AMratio)
#TODO
#Create Heyoka non-terminal event function to detect when possible collisions accur
#Create event function callback method that logs data of possible collisions when they are detected
#Write implementation of integration in PARALLEL :( (joblib library)
#
#
#Probably lots of other stuff