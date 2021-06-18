import numpy as np
from FunctionLibrary import ElementIn


list1 = np.linspace(10**-0.5, 10**1.8, 100)

#print(list1)

elementList = ElementIn("elements",".txt")

print(elementList)

