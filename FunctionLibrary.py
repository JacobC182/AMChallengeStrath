#This file contains coordinate system conversion functions
#For converting between a geocentric cartesian system and Orbital elements and opposite
from copy import Error
from numpy.lib.function_base import median
from numpy.ma.core import dot, shape
from orbital.utilities import eccentric_anomaly_from_mean, mean_anomaly_from_eccentric

#importing NumPy library
import numpy as np
from pandas.core.accessor import register_dataframe_accessor
import pykep

def vec2orb(pos, vel):
#Function for converting a state vector to orbital elements (Around Earth/Geocentric frame)
#UNITS - KM, KM/s (INPUTS)
#OUTPUT UNITS ARE: KM, RADIANS

#Defining constants
    u = np.double(3.986004407799724e+5)   #GM Earth - km3/sec2

#calculating magnitude of position vector
    mag = np.linalg.norm(pos)
#calculating orbital momentum vector
    h = np.cross(pos,vel)

#calculating eccentricity vector
    e = (np.cross(vel,h)/u) - (pos/mag)

#calculating vector that points towards ascending node
    n = [-h[1], h[0], 0]

#calculating true anomaly (radians)
    nu = np.arccos(np.dot(e,pos)/(np.linalg.norm(e)*mag))

#true anomaly quadrant check
    if np.dot(pos,vel) < 0:
        nu = 2*np.pi - nu

#calculating inclination
    i = np.arccos(h[2]/np.linalg.norm(h))

#calculating eccentricity magnitude
    emag = np.linalg.norm(e)

#calculating eccentric anomaly
    E = 2*np.arctan(np.tan(nu/2)/np.sqrt((1+emag)/(1-emag)))

#calculating longitude of the ascending node
    raan = np.arccos(n[0]/np.linalg.norm(n))

#raan quadrant check
    if n[0] < 0:
        raan = 2*np.pi - raan

#calculating argument of periapsis
    argp = np.arccos(np.dot(n,e)/(np.linalg.norm(n)*emag))

#argp quadrant check
    if e[2] < 0:
        argp = 2*np.pi - argp
    
#calculating mean anomaly
    M = E - emag*np.sin(E)

#calculating semimajor axis
    a = 1/((2/mag)-((np.linalg.norm(vel)**2)/u))

#returning orbital elements
    return a, emag, i, raan, argp, nu




def orb2vec(a, e, i, raan, argp, nu):
#Function for converting from orbital elements to state vector
#INPUTS - semimajor axis-a, eccentricity-e, inclination-i, ascending node-raan
#       - argument of periapsis-argp, mean anomaly-M
#OUTPUTS
#State vector - 2 1x3 arrays of positon and velocity respectively
#UNITS KM, Radians, KM/s 

#importing NumPy library
    import numpy as np

#defining constants
    u = np.double(3.986004407799724e+5)   #GM Earth - km3/sec2

#calculating semilatus rectum
    sl = a*(1-e**2)

#calculating magnitude (polar coordinate system)
    rm = sl/(1+e*np.cos(nu))

#calculating angular position (polar coordinate system)
    lat = argp+nu

#calculating and combining position vector
    r = [(rm*(np.cos(raan)*np.cos(lat) - np.sin(raan)*np.cos(i)*np.sin(lat))), (rm*(np.sin(raan)*np.cos(lat) - np.cos(raan)*np.cos(i)*np.sin(lat))), (rm*np.sin(i)*np.sin(lat))]

#calculating velocity vector
    vx = -np.sqrt(u/sl) * (np.cos(raan)*(e* np.sin(argp) + np.sin(lat)) + np.sin(raan)*np.cos(i)*(e*np.cos(argp)+np.cos(lat)) )
    vy = -np.sqrt(u/sl) * (np.sin(raan)*(e* np.sin(argp) + np.sin(lat)) - np.cos(raan)*np.cos(i)*(e*np.cos(argp)+np.cos(lat)) )
    vz = np.sqrt(u/sl) * (e * np.cos(argp) + np.cos(lat)) * np.sin(i)
#combining velocity vector
    v = [vx, vy, vz]

#returning position and velocity vectors
    return r, v




#callback Function for integrator
#for use as a NON-TERMINAL event with the integrator

def callback1(ta, time, dsign):
    print("callback function triggered")





#Function for taking in the debris orbital elements

def ElementIn(elementFilePath,fileExtension):

    from os.path import exists

    limiter = " "
    elements = []
    i = 0

    while exists((elementFilePath + (str(i)) + fileExtension)):

        elementLine = np.loadtxt(fname = (elementFilePath + str(i) + fileExtension), delimiter = limiter)

        elements.append(list(elementLine))

        i += 1

    return(elements)




#Collision detection callback function

def CollisionCB(ta, time, d_sign):
    print("Collision Detected")

    file = open("Collisions.txt", "a+")

    savestate = [time, ta.state[0], ta.state[1], ta.state[2], ta.state[3], ta.state[4], ta.state[5]]

    np.savetxt(file, [savestate])

    file.close()



#Function that calculates eccentric anomaly using Newton-Raphson method from the mean anomaly and eccentricity
def EccentricAnomalySolver(mean, e):
    
    import math as ma

    mean = mean * (np.pi/180)

    if e > 0.5:
        e0 = 3.14159
    else:
        e0 = mean
#while using error/corrector term condition
    while abs(e0 - (mean+e*ma.sin(e0)))/(1-e*ma.cos(e0)) > 0.0000001:
#n-r method
        e1 = e0 - ((e0 - (mean+e*ma.sin(e0)))/(1-e*ma.cos(e0)))
        print(e1/(ma.pi/180))
        e0 = e1

    return(e0)




#Function that reads debris observation files and returns the converted state vectors and time vector
def DebrisRead(fileNumber):

    from pykep import par2ic
    import numpy as np

    filename = "data\deb_train\eledebtrain" + str(fileNumber) + ".dat"

    debElements = np.loadtxt(filename)

    timeVector = []
    debVector = []

    for i in range(len(debElements)):
        eAnom = EccentricAnomalySolver(debElements[i][4], debElements[i][2])
        timeVector.append(debElements[i][0])
        debVector.append(par2ic([debElements[i][1], debElements[i][2], debElements[i][3]*(np.pi/180), debElements[i][6]*(np.pi/180), debElements[i][5]*0.0174532925199, eAnom],3.986004407799724e+5))

    return timeVector, np.reshape(debVector, [len(debElements), 6])




#Function that reads debris observation files and returns time vector and ORBITAL ELEMENTS
def DebrisReadElement(fileNumber):

    filename = "data\deb_train\eledebtrain" + str(fileNumber) + ".dat"

    debElements = np.loadtxt(filename)

    timeVector = debElements[:,0]
    debVector = debElements[:, 1:-1:1]

    return timeVector, debVector




#Function that converts state vector to orbital elements - INPUT = 1x6 state vector
def rv2orb(state):

    from pykep import ic2par

    r = state[0:3]

    v = state[3:6]

    orb = ic2par(r,v,3.986004407799724e+5)

    return orb

#Function that converts state vector to orbital elements - INPUT = 1x6 state vector - OUTPUTS in observation file format
def rv2orbF(state):

    from pykep import ic2par

    r = state[0:3]

    v = state[3:6]

    orb = ic2par(r,v,3.986004407799724e+5)

    MeanAnom = orb[5]  - (orb[1]*np.sin(orb[5]))
    #MeanAnom = mean_anomaly_from_eccentric(orb[1],orb[5])

    orb = [orb[0], orb[1], orb[2]*(180/np.pi), (MeanAnom*(180/np.pi))%360, orb[4]*(180/np.pi), orb[3]*(180/np.pi)]

    return orb




#Function that reads and returns the initial state information from the labels-training file
def DebrisLabel(fileNumber):

    n = int(fileNumber) - 1 #convert filenumber specifier from string to integer (Why did you have to list the files with leading zeros!?!)

    bigList = np.loadtxt("data\labels_train.dat")   #reading entire labels_train file

    return bigList[n,:] #returning only the row corresponding to the chosen debris file number




#Function that creates and trains the AREA_TO_MASS ratio Machine Learning Model! - Uses the Random Forest-Decision Tree Regressor Model (Not Anymore) ->Gradient Boosting Regressor
def CRAMRegressorModel():

    #Importing Gradient Boost Regressor (Black Box) Class - Used to instantiate Regressor Model Objects - From Scikit-Learn (SKlearn) Machine Learning Library
    from sklearn.ensemble import GradientBoostingRegressor
    #Creating Gradient Boosting Regressor (Black Box) Object
    Regressor = GradientBoostingRegressor(n_estimators=2000, max_depth=2)

    X = []      #Create empty dataset arrays
    AMratio = []

    #Reading training data output values into file
    debrisTrainRatio = np.loadtxt("data\labels_train.dat")[:,1]
    DebrisTrainState = np.loadtxt("data\labels_train.dat")[:,2:9]
    
#READING DATA----------------------------------------------------------------------------------
#Iterating over a range 1-100 step=1, going through all 100 debris data files
    for i in range(26, 100 +1, 1):
        fileNumberString = ""   #creating empty string for converting int (i) to string to use as sequential file number

        if len(str(i)) == 1:    #Creating appropriate string of file number "001" - "100"
            fileNumberString = "00" + (str(i))
        elif len(str(i)) == 2:
            fileNumberString = "0" + (str(i))
        else:
            fileNumberString = (str(i))
    
        debrisData = np.loadtxt("data\deb_train\eledebtrain" + fileNumberString + ".dat")   #Reading each debris observation file individually

        if len(np.shape(debrisData)) == 1:      #Appending debris data and AM ratio to input and output training arrays respectivel
            AMratio.append(debrisTrainRatio[i-1])
            X.append(debrisData)
        else:
            for j in range(len(debrisData)):        #Multiple observations from the same debris are given the same true AM-ratio from that debris in the 2D training array format
                AMratio.append(debrisTrainRatio[i-1])
                X.append(debrisData[j,:])
#READING DATA END------------------------------------------------------------------------------

    Regressor.fit(X, AMratio)   #Training "Black Box" Regressor to training data

    print("Fitting Score: " + str(Regressor.score(X, AMratio)) )
    #print("Feature Importances")
    #print(Regressor.feature_importances_)
    
    return Regressor    #Returning Trained Regressor Model Object





#Function that reads ALL satellite trajectory data and returns it as a 3D array, First axis/index = satellite number
def SatRead():

    allSatData = []

    for i in range(1, 100 +1, 1):
        fileNumberString = ""   #creating empty string for converting int (i) to string to use as sequential file number

        if len(str(i)) == 1:    #Creating appropriate string of file number "001" - "100"
            fileNumberString = "00" + (str(i))
        elif len(str(i)) == 2:
            fileNumberString = "0" + (str(i))
        else:
            fileNumberString = (str(i))
    
        SatData = np.loadtxt("data\sat\elesat" + fileNumberString + ".dat")     #Reading individual satellite trajectory file
        allSatData.append(SatData)      #Adding SatData (read-in file) to 3D array of all trajectory data

    return allSatData       #Returning complete 3d array of all satellite data





#Function that returns the string of the input integer from 1-100, with a str-length of 3, filled with leading zeros
#eg 1 = "001", 57 = "057"
def FileStr(i):

    fileNumberString = ""   #creating empty string for converting int (i) to string to use as sequential file number

    if len(str(i)) == 1:    #Creating appropriate string of file number "001" - "100"
        fileNumberString = "00" + (str(i))
    elif len(str(i)) == 2:
        fileNumberString = "0" + (str(i))
    else:
        fileNumberString = (str(i))
    
    return fileNumberString




#Function that returns the MEAN SQUARE ERROR from a set of "actual" and "predicted" values - TAKEN FROM: https://www.statology.org/mean-squared-error-python/
def mse(actual, pred): 
    actual, pred = np.array(actual), np.array(pred)
    return np.square(np.subtract(actual,pred)).mean()




#Function that converts system of orbital elements taken from data files in form a,e,i,M,w,W, returns 1x0 state vector in KM, KM/s, around Earth
def orb2rv(startState):

    #startR, startV = pykep.par2ic([startState[0], startState[1], startState[2]*(np.pi/180), startState[5]*(np.pi/180), startState[4]*(np.pi/180), EccentricAnomalySolver(startState[3], startState[1])], 3.986004407799724e+5)
    startR, startV = pykep.par2ic([startState[0], startState[1], startState[2]*(np.pi/180), startState[5]*(np.pi/180), startState[4]*(np.pi/180), eccentric_anomaly_from_mean(startState[1], startState[3]*(np.pi/180))], 3.986004407799724e+5)

    return np.concatenate((startR, startV), axis=0)






#Function that creates and trains the Observational Outlier Detector Model!! - sklearn outlier detection - machine learning model
def XrayVision():
    
    #importing outlier detection algorithm models from sklearn library
    from sklearn.ensemble import IsolationForest
    from sklearn.neighbors import LocalOutlierFactor
    #Unsure if neighbors algorithm is useable - "Curse of dimensionality, Concept of neighborhood becomes meaningless" --https://archive.siam.org/meetings/sdm10/tutorial3.pdf - Slide 56
    #But neighbors model accepts N x M array multi-dimension data
    from sklearn.svm import OneClassSVM
    from sklearn.covariance import EllipticEnvelope

    #Creating Outlier Detection Model object
    OutlierModel = IsolationForest(n_estimators=2000, n_jobs=-1)
    #OutlierModel = LocalOutlierFactor(n_neighbors=2)
    
    #Creating empty array to hold debris data
    debData = []

    #READING DEBRIS DATA-------------------------------------------------
    for i in range(26, 100 +1, 1):
        fileNumberString = ""   #creating empty string for converting int (i) to string to use as sequential file number

        if len(str(i)) == 1:    #Creating appropriate string of file number "001" - "100"
            fileNumberString = "00" + (str(i))
        elif len(str(i)) == 2:
            fileNumberString = "0" + (str(i))
        else:
            fileNumberString = (str(i))
    
        debrisData = np.loadtxt("data\deb_train\eledebtrain" + fileNumberString + ".dat")   #Reading each debris observation file individually

        if len(np.shape(debrisData)) == 1:      #Appending debris data and AM ratio to input and output training arrays respectivel
            debData.append(debrisData)
        else:
            for j in range(len(debrisData)):        #Multiple observations from the same debris are given the same true AM-ratio from that debris in the 2D training array format
                debData.append(debrisData[j,:])

    #Fitting Model to training data
    OutlierModel.fit(debData)

    return OutlierModel






#XrayVision()