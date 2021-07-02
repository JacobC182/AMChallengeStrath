import numpy as np
import heyoka as hy
from starterODE import ODE
import scipy.optimize
sysODE = ODE()

t0 = -1935.00*60*60*24
initial = [-3.84481241e+04, -1.13678964e+04, -1.76824887e+03,  9.17589734e-01, -3.00803552e+00, -7.38669271e-01]
t1 = -255.00*60*60*24
final = [3.17895790e+04,  1.33646409e+04,  3.47242170e+03, -1.43104345e+00, 3.30468682e+00,  7.77195904e-01]
AM0 = 28*1e-6

ta = hy.taylor_adaptive(sys = sysODE, state = initial)

Amlimit = scipy.optimize.Bounds((10**-0.5)*1e-6, (10**1.8)*1e-6)


def ObjFunction(AM):

    ta.state[:] = initial
    ta.time = t0
    ta.pars[0] = AM

    ta.propagate_until(t=t1)

    out = ta.state

    min = np.subtract(final[0:3], out[0:3])
    print(out[0:3], final[0:3])
    return(np.linalg.norm(min))

optimal = scipy.optimize.minimize(ObjFunction, [AM0], bounds = Amlimit, method = "TNC")

print(optimal)

