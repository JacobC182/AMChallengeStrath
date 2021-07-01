from os import write
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from FunctionLibrary import *

#Creating Random Forest Regressor (Black Box) Object
Regressor = RandomForestRegressor()


#Creating empty X and Y arrays
X = []
AMratio = []

#Reading training data output values into file
debrisTrainRatio = np.loadtxt("data\labels_train.dat")[:,1]



for i in range(1, 101, 1):
    

    fileNumberString = ""

    if len(str(i)) == 1:
        fileNumberString = "00" + (str(i))
    elif len(str(i)) == 2:
        fileNumberString = "0" + (str(i))
    else:
        fileNumberString = (str(i))

    debrisData = np.loadtxt("data\deb_train\eledebtrain" + fileNumberString + ".dat")

    if len(np.shape(debrisData)) == 1:

        AMratio.append(debrisTrainRatio[i-1])

        X.append(debrisData)
    else:

        for j in range(len(debrisData)):
        
            AMratio.append(debrisTrainRatio[i-1])

            X.append(debrisData[j,:])
    
    print("Finished Reading File: " + fileNumberString)

print("Starting Training")
Regressor.fit(X, AMratio)

#wronganswer = Regressor.predict([[-1935.00,   42276.74,  0.0507489,    13.4614,   107.7131,   188.8823, 5.8539]])
wronganswer = Regressor.predict([[-235.00  , 42249.88 , 0.1497848  ,  20.2557  ,  79.3143  , 360.5050  , 356.3467]])



print(wronganswer)