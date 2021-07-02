import numpy as np


lengths = []
lengthsTrain = []
for i in range(100):
    fileNum = i + 1

    if len(str(fileNum)) == 1:
        fileNumStr = "00" + str(fileNum)
    elif len(str(fileNum)) == 2:
        fileNumStr = "0" + str(fileNum)

    debFile = np.loadtxt("data\deb_test\eledebnewfd" + fileNumStr + ".dat")
    debFileTrain = np.loadtxt("data\deb_train\eledebtrain" + fileNumStr + ".dat")

    debFile = np.reshape(debFile, [-1,1])
    debFileTrain = np.reshape(debFile, [-1,1])

    lengths.append(len(debFile)/7)
    lengthsTrain.append(len(debFileTrain)/7)

only1 = 0
only1Train = 0

for i in lengths:
    if i == 1:
        only1 +=1

for i in lengthsTrain:
    if i == 1:
        only1Train +=1

print("No. of Debris with only 1 observation: " + str(only1) )
print("No. of Training Debris with only 1 observation: " + str(only1Train) )
print("---")
print("Total Debris Observations: " + str(sum(lengths)))
print("Total Training Debris Observations: " + str(sum(lengthsTrain)))