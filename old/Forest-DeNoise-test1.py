from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
import numpy as np

ErrorRegressor = RandomForestRegressor(n_estimators=1000, n_jobs=-1)

observedInput = np.loadtxt("observedData.txt",delimiter=",")
propagatedOutput = np.loadtxt("propagatedData.txt",delimiter=",")

ErrorRegressor.fit(observedInput, propagatedOutput)

deNoise = ErrorRegressor.predict([[42276.74 , 0.0507489  ,  13.4614  , 107.7131 ,  188.8823  ,   5.8539]])

print(deNoise)