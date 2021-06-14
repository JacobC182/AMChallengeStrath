from Integrator import int1
from ODE import ODE

coreCount = 8

initial = [2.6533e4, 0, 0, 0, 2.2220, 3.1734]

ODE = ODE()

import multiprocessing as mp

#p1 = mp.Process(target=int1(ODE, initial, 200000, [0.1]))
#p2 = mp.Process(target=int1(ODE, initial, 200000, [0.2]))
#p3 = mp.Process(target=int1(ODE, initial, 200000, [0.3]))
#p4 = mp.Process(target=int1(ODE, initial, 200000, [0.4]))
#p5 = mp.Process(target=int1(ODE, initial, 200000, [0.5]))
#p6 = mp.Process(target=int1(ODE, initial, 200000, [0.6]))
#p7 = mp.Process(target=int1(ODE, initial, 200000, [0.7]))
#p8 = mp.Process(target=int1(ODE, initial, 200000, [0.8]))

AM = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]

from multiprocessing_on_dill import Pool

if __name__ == "__main__":
    pool1=Pool(processes=3)

    r = pool1.map(int1,[ODE,initial,200000,0.1] )

#THIS DOES NOT WORK
#ParallelPython or multiprocessing, or joblib cannot successfully interact with the Heyoka library - "heyoka.core.expression" objects