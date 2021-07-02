import time
import heyoka as hy
from numba.core.types.misc import EllipsisType
import numpy as np
from FunctionLibrary import *
from starterODE import ODE
import pykep as pk
import matplotlib.pyplot as plt

t, v = DebrisRead("001")

initial = v[0]

ta = hy.taylor_adaptive(sys = ODE(), state = initial, pars = [10.32583*1e-6])

tgrid = np.linspace(start = (t[0] * pk.SEC2DAY), stop = (t[0] * pk.SEC2DAY)+365.25*24*60*60*5, num = 60*60*24, endpoint=True)

out = ta.propagate_grid(tgrid)

vec = out[4]


GMe = 3.986004407799724e+5
size = np.shape(vec)
print(size)
ele = np.zeros(size)

for i, cart in enumerate(vec):
    ele[i] = pk.ic2par(cart[:3], cart[3:], GMe)


print(ele)

fig, axs = plt.subplots(3, 2)

axs = axs.reshape(6,)
for i in range(6):
    axs[i].plot(tgrid * pk.SEC2DAY, ele[:,i]/ele[:,2])
axs[0].set_ylabel("a (km)")
axs[1].set_ylabel("e")
axs[2].set_ylabel("i (rad)")
axs[3].set_ylabel("W (rad)")
axs[4].set_ylabel("w (rad)")
axs[5].set_ylabel("E (rad)")

plt.tight_layout()

plt.show()

