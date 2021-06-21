import time
import heyoka as hy
import numpy as np
import math as ma
import pykep as pk
from joblib.parallel import delayed
from starterODE import ODE
from FunctionLibrary import *
import scipy.optimize as optimization

#Starting script timer

#Creating ODE system
sysODE = ODE()

#Reading debris data file
debElements = np.loadtxt("data\deb_train\eledebtrain001.dat")

debVector = []

for i in range(len(debElements)):

    eAnom = EccentricAnomalySolver(debElements[i][4], debElements[i][2])

    debVector.append(pk.par2ic([debElements[i][1], debElements[i][2], debElements[i][3], debElements[i][6], debElements[i][5], eAnom]))

