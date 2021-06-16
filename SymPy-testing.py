#This is a TEST - that failed
#This script aided in the tests where SymPy was used to simplify the Kep, J2, C22, and S22 equations

#simplifying ODE system with SymPy module
import heyoka as hy
from ODE import ODE
import sympy as sm

ODE1 = ODE([1,1,1,1,1,1,1])
#print(ODE1[3][1])
smODE = hy.to_sympy(ODE1[5][1])

simpleODE_sympy = sm.simplify(smODE)

#print(simpleODE_sympy)
#print(hy.from_sympy(simpleODE_sympy))

#saving simplified "sympy" equation to file
e1 = open('e3.txt', "a")
e1.write(str(simpleODE_sympy))
e1.close