from joblib.parallel import delayed
from numpy import core
from ODE import ODE
import numpy as np
import time
import heyoka as hy

#setting time datum
startTime = time.time()

coreCount = 8

initial = [2.6533e4, 0, 0, 0, 2.2220, 3.1734]

ODE = ODE()

AM = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]

import joblib as jl
#joblib appears to work when setting the preference to "threads"

ta = hy.taylor_adaptive(sys=ODE, state=initial)

def int1(sysODE, init,t, AM):

    ta.pars[0] = AM
    ta.time = 0
    ta.state[:] = init

    ta.propagate_until(t)

    return ta.state

output = jl.Parallel(n_jobs=coreCount,prefer="threads")(delayed(int1)(ODE,initial,2000,AM[i]) for i in range(8))
print(output)
#np.savetxt("parOut.txt",output,delimiter=",")

print("Script Run Time: " + str(time.time() - startTime)[0:8] + "s")