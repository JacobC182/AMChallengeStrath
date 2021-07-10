#Script is used to test a multitude of different machine learning algorithms and conf9igurations to try and get an optimal setup
from os import write
from sklearn.svm import SVR
import numpy as np
from sklearn.model_selection import RandomizedSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import SGDRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import GradientBoostingRegressor
#Creating Random Forest Regressor (Black Box) Object
#Regressor = SVR(kernel="poly", degree=80, C=0.9, tol=1e-5, cache_size=4096)
#Regressor = DecisionTreeRegressor()
#Regressor = KNeighborsRegressor(n_jobs=-1, n_neighbors=2, leaf_size=30, )
#Regressor = KernelRidge(alpha=0.9)
#Regressor = MLPRegressor(hidden_layer_sizes=(600,), solver="adam", learning_rate="adaptive", max_iter=350)
Regressor = GradientBoostingRegressor(n_estimators=3000, tol=1e-8)

#Creating empty X and Y arrays
X = []
AMratio = []

#Reading training data output values into file
debrisTrainRatio = np.loadtxt("data\labels_train.dat")[:,1]

print("Reading Training Dataset")
#Iterating over a range 1-100 step=1, going through all 100 debris data files
for i in range(1, 98 +1, 1):
    

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

sample = np.loadtxt("data\deb_train\eledebtrain099.dat")
sample = np.reshape(sample, [-1,7])

print("Finished Training, Starting Prediction")
wronganswer = Regressor.predict(sample)       #Predicting output value from input sample data using trained model
print(wronganswer)
wronganswer = sum(wronganswer)/len(wronganswer)

print("Prediction:  " + str(wronganswer) )
#print(Regressor._get_param_names)
print(Regressor.score(X, AMratio))

#TUNING------

#rTune = RandomizedSearchCV(estimator=Regressor, )