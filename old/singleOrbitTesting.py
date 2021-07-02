import heyoka as hy
import numpy as np
from FunctionLibrary import *
from starterODE import ODE

sysODE = ODE()

initial = [-3.84481241e+04, -1.13678964e+04, -1.76824887e+03, 9.17589734e-01, -3.00803552e+00, -7.38669271e-01]

AMratio = 10.32583*1e-6

ta = hy.taylor_adaptive(sys = sysODE, state = initial, pars = [AMratio], high_accuracy=True, tol = 1e-17)

ta.time = -1935.00*60*24*60


status, min_h, max_h, nsteps = ta.propagate_until(t=-255.00*24*60*60)

print("Outcome      : {}".format(status))
print("Min. timestep: {}".format(min_h))
print("Max. timestep: {}".format(max_h))
print("Num. of steps: {}".format(nsteps))
print("Current time : {}\n".format(ta.time))
print(ta)

correct = [3.17895790e+04,  1.33646409e+04, 3.47242170e+03, -1.43104345e+00, 3.30468682e+00, 7.77195904e-01]