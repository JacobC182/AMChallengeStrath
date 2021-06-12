from Integrator import int1
from ODE import ODE

coreCount = 8

initial = [2.6533e4, 0, 0, 0, 2.2220, 3.1734]

import pp

job_server = pp.Server(ncpus=coreCount)
print ("Starting Parallel Pool with " + str(coreCount) + "Cores")


for i in range(coreCount):
    job_server.submit(int1, args=(ODE(),initial,200000,[0.1]), modules=("heyoka","numpy","math"))

#THIS DOES NOT WORK
#ParallelPython cannot successfully interact with the Heyoka library - "heyoka.core.expression" objects