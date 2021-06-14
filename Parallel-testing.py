from joblib.parallel import delayed
from numpy import core
from Integrator import int1
from ODE import ODE
import numpy as np
import time

#setting time datum
startTime = time.time()

coreCount = 8

initial = [2.6533e4, 0, 0, 0, 2.2220, 3.1734]

ODE = ODE()

AM = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]

import joblib as jl
#joblib appears to work when setting the preference to "threads"

output = jl.Parallel(n_jobs=coreCount,prefer="threads")(delayed(int1)(ODE,initial,20000,[AM[i]]) for i in range(8))

np.savetxt("parOut.txt",output,delimiter=",")

print("Script Run Time: " + str(time.time() - startTime)[0:8] + "s")