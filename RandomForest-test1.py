from os import write
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from FunctionLibrary import *

#Creating Random Forest Regressor (Black Box) Object
Regressor = RandomForestRegressor(n_estimators=2000, criterion="mse", n_jobs=-1, verbose=1)


#Creating empty X and Y arrays
X = []
AMratio = []

#Reading training data output values into file
debrisTrainRatio = np.loadtxt("data\labels_train.dat")[:,1]

print("Reading Training Dataset")
#Iterating over a range 1-100 step=1, going through all 100 debris data files
for i in range(1, 100 +1, 1):
    

    fileNumberString = ""   #creating empty string for converting int (i) to string to use as sequential file number

    if len(str(i)) == 1:    #Creating appropriate string of file number "001" - "100"
        fileNumberString = "00" + (str(i))
    elif len(str(i)) == 2:
        fileNumberString = "0" + (str(i))
    else:
        fileNumberString = (str(i))

    debrisData = np.loadtxt("data\deb_train\eledebtrain" + fileNumberString + ".dat")   #Reading each debris observation file individually

    if len(np.shape(debrisData)) == 1:      #Appending debris data and AM ratio to input and output training arrays respectively

        AMratio.append(debrisTrainRatio[i-1])

        X.append(debrisData)
    else:

        for j in range(len(debrisData)):        #Multiple observations from the same debris are given the same true AM-ratio from that debris in the 2D training array format
        
            AMratio.append(debrisTrainRatio[i-1])

            X.append(debrisData[j,:])
    

print("Starting Training")
Regressor.fit(X, AMratio)   #Training "Black Box" Regressor to training data

sample = np.loadtxt("data\deb_train\eledebtrain100.dat")

print("Finished Training, Starting Prediction")
wronganswer = Regressor.predict([sample])       #Predicting output value from input sample data using trained model


print("Prediction:  " + str(wronganswer) )