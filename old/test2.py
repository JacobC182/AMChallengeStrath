import numpy as np
from FunctionLibrary import *

timeArray = []

for i in range(1, 100 +1, 1):
    
    
    fileNumberString = ""   #creating empty string for converting int (i) to string to use as sequential file number

    if len(str(i)) == 1:    #Creating appropriate string of file number "001" - "100"
        fileNumberString = "00" + (str(i))
    elif len(str(i)) == 2:
        fileNumberString = "0" + (str(i))
    else:
        fileNumberString = (str(i))

    debrisData = np.loadtxt("data\deb_train\eledebtrain" + fileNumberString + ".dat")
    debrisData = np.reshape(debrisData, [-1,7])
    timeArray = np.append(timeArray, debrisData[:,0])

    timeArray = np.reshape(timeArray, [-1,1])

bigData = np.loadtxt("propagatedData.txt", delimiter = ",")

biggerData = np.concatenate((timeArray, bigData), axis = 1)

np.savetxt(fname="PropData", X=biggerData)