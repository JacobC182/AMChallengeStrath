import numpy as np
import heyoka as hy
from ODE import ODE
from FunctionLibrary import orb2rv, rv2orbF

sysODE = ODE()

ta = hy.taylor_adaptive(sys=sysODE, state = [0, 0, 0, 0, 0, 0])

def Int(startTime, endTime, startState, AMratio):

    ta.time = startTime *60*60*24
    ta.state[:] = startState
    ta.pars[0] = AMratio *1e-6

    ta.propagate_until(endTime *60*60*24)

    return ta.state[:]

choice = 2

AMsample = [10.32583, 20.08288, 0.52148, 12.81137]
debSample = [[-1935, 42276.74, 0.0507489, 13.4614, 107.7131, 188.8823, 5.8539], [-235.00, 42271.82, 0.0637583, 27.2531, 112.5560, 338.9246, 348.3112], [3695.00, 42284.33, 0.1011880, 16.6487, 256.3106, 261.1312, 348.9994], [-235.00, 42249.88, 0.1497848, 20.2557, 79.3143, 360.5050, 356.3467]]
endtimeSample = [-3815, -3845, -7275, -3485]

statevec = orb2rv(debSample[choice][1:7])


outstate = Int(startTime=debSample[choice][0], endTime=endtimeSample[choice], startState=statevec, AMratio=AMsample[choice])

#print(outstate)
print(rv2orbF(outstate))
print("End Time: " + str(ta.time /(60*60*24)) )

"""     1,2,3,4 - sat 94,17,18,70
42260.01  0.0477003    11.3705    25.6240    91.5576     9.5178
42260.66  0.0081106     8.2732    56.3063   272.0286    45.9370
42252.96  0.0996132     2.2224    13.8162   242.1460   166.9770
42253.07  0.0077042    14.8392    54.9740   203.1057    30.9312
"""