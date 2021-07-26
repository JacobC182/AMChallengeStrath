import numpy as np

residual = np.loadtxt("residual.txt")

outliers = 0

for i in range(len(residual)):
    if abs(residual[i,3]) > 10 or abs(residual[i,4]) > 10 or abs(residual[i,5]) > 10:
        outliers += 1

print("No. of outliers: " + str(outliers))