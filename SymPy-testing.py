#simplifying ODE system with SymPy module
import heyoka as hy
from ODE import ODE
import sympy
ODE1 = ODE([1,1,1,1,1,1,1])

smODE = hy.to_sympy(ODE1)


def SymPySimplify(smODE):
    import sympy as sm
    return hy.from_sympy(sm.simplify(smODE))

simpleODE = SymPySimplify(smODE)

print(simpleODE)