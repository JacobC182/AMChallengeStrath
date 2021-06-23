import heyoka as hy
import numpy as np
from FunctionLibrary import *
from starterODE import ODE

sysODE = ODE()

initial = [40538.01113051619, 10803.588965445128, -3143.1486811323193, -0.7408330964852461, 2.9707963804406887, 0.37782171710291634]

AMratio = 10.32583*1e-3

ta = hy.taylor_adaptive(sys = sysODE, state = initial, pars = [AMratio])

ta.time = -3815.00*60*24*60


status, min_h, max_h, nsteps = ta.propagate_until(t=-1935.00*24*60*60)

print("Outcome      : {}".format(status))
print("Min. timestep: {}".format(min_h))
print("Max. timestep: {}".format(max_h))
print("Num. of steps: {}".format(nsteps))
print("Current time : {}\n".format(ta.time))
print(ta)