import numpy as np

from FunctionLibrary import FileStr, mse
from sklearn.ensemble import IsolationForest, GradientBoostingRegressor
from sklearn.svm import OneClassSVM
from sklearn.covariance import EllipticEnvelope

data = np.loadtxt("DebDataCompiled.txt")
AMratio = np.loadtxt("AMRatioCompiled.txt")
debTrain = data[20:-1, :]
debTest = data[0:20 +1, :]


debTrain = np.reshape(debTrain, [-1, 7])
debTest = np.reshape(debTest, [-1, 7])


FogLight = IsolationForest(n_estimators=200, n_jobs=-1)#, max_samples=200)
FogLight = OneClassSVM(nu=0.3)
#FogLight = EllipticEnvelope()

FogLight.fit(debTrain)

results = FogLight.fit_predict(debTest)
print(FogLight.support_vectors_)
print(results)

Regressor = GradientBoostingRegressor(n_estimators=2000, max_depth=2)
Regressor.fit(debTrain, AMratio[20:-1])
print("Fitting Score: " + str(Regressor.score(debTrain, AMratio[20:-1])) )

AMratioOut = []
for i in range(0,20, 1):
    predictedRatio = Regressor.predict(debTest)     #Predicting AM-ratio using Model

    predictedRatio = sum(predictedRatio)/len(predictedRatio)    #Taking average of predictions for multiple observation debris
    AMratioOut.append(predictedRatio)

realRatio = np.loadtxt("data\labels_train.dat")[:,1]
msError = mse(realRatio[0:20], AMratioOut)
print("MSE: " + str(msError))