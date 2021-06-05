#2021 Andrea Milani Challenge
#THIS IS A TEST GITHUB IS SHITE
#Importing time module - to measure script execution time
import time
scriptStartTime = time.time()

#Importing "Heyoka" integrator library & NumPy library & Math library
import heyoka as hy
import numpy as np
import math as ma

#Defining Area-to-Mass Ratio (TEMPORARY PLACEHOLDER)
Aratio = 1

#Defining Constants
GME = np.double(3.986004407799724e+5)   #GM Earth - km3/sec2
GMS = np.double(1.32712440018e+11)      #GM Sun - km3/sec2
GMM = np.double(4.9028e+3)              #GM Moon - km3/sec2
RE = np.single(6378.1363)               #Earth radius - km
C20 = np.double(-4.84165371736e-4)      #C20 coefficient
C22 = np.double(2.43914352398e-6)       #C22 coefficient
S22 = np.double(-1.40016683654e-6)      #S22 coefficient
thetaG = np.double((np.pi/180)*280.4606)
vE = np.double((np.pi/180)*280.4606)
vS = np.double((np.pi/180)*1.1407410259335311e-5)
vMa = np.double((np.pi/180)*1.512151961904581e-4)
vMp = np.double((np.pi/180)*1.2893925235125941e-6)
vMs = np.double((np.pi/180)*6.128913003523574e-7)
aS = np.double(1.49619e+8)
eps = np.double((np.pi/180)*23.4392911)
phiS0 = np.double((np.pi/180)*357.5256)
OmegaS = np.double((np.pi/180)*282.94)
pSRP = np.double(4.56e-6)

#Defining shortcut variables (Ease of Use)
r15 = ma.sqrt(15)                       #Square root of 15
r5 = ma.sqrt(5)                         #Square root of 5
t = hy.time                             # "t" defined as the time of the integrator "heyoka.time"

#Creating symbolic variables
X, Y, Z = hy.make_vars("X","Y","Z")
pX, pY, pZ = hy.make_vars("pX","pY","pZ")
#Defining Equations of motion

#x,y,z terms for C22 and S22 terms
x = X*hy.cos(thetaG+(vE*t)) + Y*hy.sin(thetaG+(vE*t))
y = -X*hy.sin(thetaG+(vE*t)) + Y*hy.cos(thetaG+(vE*t))
z = Z

#C22 terms
fC22x = ( (5*GME*(RE**2)*r15*C22*x*((y**2)-(x**2))) /  (2*((x**2+y**2+z**2)**3.5)) ) + ( (GME*(RE**2)*r15*C22*x)/((x**2+y**2+z**2)**2.5) )
fC22y = ( (5*GME*(RE**2)*r15*C22*y*((y**2)-(x**2))) /  (2*((x**2+y**2+z**2)**3.5)) ) - ( (GME*(RE**2)*r15*C22*y)/((x**2+y**2+z**2)**2.5) )
fC22z = (5*GME*(RE**2)*r15*C22*z*((y**2)-(x**2))) /  (2*((x**2+y**2+z**2)**3.5))

#S22 terms
fS22x = ( (-5*GME*(RE**2)*r15*S22*(x**2)*y) /  ((x**2+y**2+z**2)**3.5) ) + ( (GME*(RE**2)*r15*S22*y)/((x**2+y**2+z**2)**2.5) )
fS22y = ( (-5*GME*(RE**2)*r15*S22*(y**2)*x) /  ((x**2+y**2+z**2)**3.5) ) + ( (GME*(RE**2)*r15*S22*x)/((x**2+y**2+z**2)**2.5) )
fS22z = ( (-5*GME*(RE**2)*r15*S22*x*y*z) /  ((x**2+y**2+z**2)**3.5) )

#C22 final terms
fC22X = fC22x*hy.cos(thetaG+(vE*t)) - fC22y*hy.sin(thetaG+(vE*t))
fC22Y = fC22x*hy.sin(thetaG+(vE*t)) + fC22y*hy.cos(thetaG+(vE*t))
fC22Z = fC22z

#S22 final terms
fS22X = fS22x*hy.cos(thetaG+(vE*t)) - fS22y*hy.sin(thetaG+(vE*t))
fS22Y = fS22x*hy.sin(thetaG+(vE*t)) + fS22y*hy.cos(thetaG+(vE*t))
fS22Z = fS22z

#J2 terms
fJ2X = ( (GME*(RE**2)*r5*C20*X)/(2*((X**2+Y**2+Z**2)**0.5)) )* ( (3/((X**2+Y**2+Z**2)**2)) - (15*Z**2)/((X**2+Y**2+Z**2)**3) )
fJ2Y = ( (GME*(RE**2)*r5*C20*Y)/(2*((X**2+Y**2+Z**2)**0.5)) )* ( (3/((X**2+Y**2+Z**2)**2)) - (15*Z**2)/((X**2+Y**2+Z**2)**3) )
fJ2Z = ( (GME*(RE**2)*r5*C20*Z)/(2*((X**2+Y**2+Z**2)**0.5)) )* ( (9/((X**2+Y**2+Z**2)**2)) - (15*Z**2)/((X**2+Y**2+Z**2)**3) )

#Keplerian "Kep" terms
fKepX = (-GME*X)/((X**2+Y**2+Z**2)**1.5)
fKepY = (-GME*Y)/((X**2+Y**2+Z**2)**1.5)
fKepZ = (-GME*Z)/((X**2+Y**2+Z**2)**1.5)

#Solar tide - sun angular terms
lSun = phiS0 + vS*t
rSun = 149.619 - 2.499*hy.cos(lSun) - 0.021*hy.cos(2*lSun)  #10^6KM !!
lambdaSun = OmegaS + lSun + (ma.pi/180)*((6892/3600)*hy.sin(lSun) + (72/3600)*hy.sin(2*lSun))

#Solar tide - Sun XYZ position terms
Xs = rSun*hy.cos(lambdaSun)
Ys = rSun*hy.sin(lambdaSun)*ma.cos(eps)
Zs = rSun*hy.sin(lambdaSun)*ma.sin(eps)

#Solar tide final terms
fSunX = -GMS*( ( (X-Xs)/ ( ((X-Xs)**2 +(Y-Ys)**2 +(Z-Zs)**2)**1.5) ) +  (Xs/((Xs**2+Ys**2+Zs**2)**1.5)) )
fSunY = -GMS*( ( (Y-Ys)/ ( ((X-Xs)**2 +(Y-Ys)**2 +(Z-Zs)**2)**1.5) ) +  (Ys/((Xs**2+Ys**2+Zs**2)**1.5)) )
fSunZ = -GMS*( ( (Z-Zs)/ ( ((X-Xs)**2 +(Y-Ys)**2 +(Z-Zs)**2)**1.5) ) +  (Zs/((Xs**2+Ys**2+Zs**2)**1.5)) )

#Lunar tide - moon time dependant angular terms
phiM = vS*t
phiMa = vMa*t
phiMp = vMp*t
phiMS = vMs*t
L0 = phiMp + phiMa + (ma.pi/180)*218.31617
lM = phiMa + (ma.pi/180)*134.96292
l2M = lSun + phiM + (ma.pi/180)*357.5256    #unsure if error this equation!
FM = phiMp + phiMa + phiMS + (ma.pi/180)*93.27283
DM = phiMp + phiMa + phiM + (ma.pi/180)*297.85027

#Lunar tide - moon angular terms
rMoon = 385000 -20905*hy.cos(lM) -3699*hy.cos(2*DM-lM) -2956*hy.cos(2*DM) -570*hy.cos(2*lM) +246*hy.cos(2*lM-2*DM) -205*hy.cos(l2M-2*DM) -171*hy.cos(lM+2*DM) -152*hy.cos(lM+l2M-2*DM)
lambdaMoon = L0 + (ma.pi/180)*((22640/3600)*hy.sin(lM) +(769/3600)*hy.sin(2*lM) -(4856/3600)*hy.sin(lM-2*DM) +(2370/3600)*hy.sin(2*DM) -(668/3600)*hy.sin(l2M) -(412/3600)*hy.sin(2*FM) -(212/3600)*hy.sin(2*lM-2*DM) -(206/3600)*hy.sin(lM+l2M-2*DM) + (192/3600)*hy.sin(lM+2*DM) -(165/3600)*hy.sin(lM-2*DM) +(148/3600)*hy.sin(lM+l2M) -(125/3600)*hy.sin(DM) -(110/3600)*hy.sin(lM+l2M) -(55/3600)*hy.sin(2*FM-2*DM))
betaMoon = (ma.pi/180)*( (18520/3600)*hy.sin(FM+lambdaMoon-L0+ ( (ma.pi/180)*( (412/3600)*hy.sin(2*FM)+(541/3600)*hy.sin(l2M))))) - (526/3600)*hy.sin(FM-2*DM) + (44/3600)*hy.sin(lM+FM-2*DM) - (31/3600)*hy.sin(-lM+FM-2*DM) - (25/3600)*hy.sin(-2*lM+FM) - (23/3600)*hy.sin(l2M+FM-2*DM) + (21/3600)*hy.sin(-lM+FM) + (11/3600)*hy.sin(-l2M+FM-2*DM)

#Lunar tide - moon XYZ terms
Xm = rMoon*hy.cos(lambdaMoon)*hy.cos(betaMoon)
Ym = (ma.cos(eps)*rMoon*hy.sin(lambdaMoon)*hy.cos(betaMoon)) + (-ma.sin(eps)*rMoon*hy.sin(betaMoon))
Zm = (ma.sin(eps)*rMoon*hy.sin(lambdaMoon)*hy.cos(betaMoon)) + (ma.cos(eps)*rMoon*hy.sin(betaMoon))

#Lunar tide final terms
fMoonX = -GMM*( ( (X-Xm)/ ( ((X-Xm)**2 +(Y-Ym)**2 +(Z-Zm)**2)**1.5) ) +  (Xm/((Xm**2+Ym**2+Zm**2)**1.5)) )
fMoonY = -GMM*( ( (Y-Ym)/ ( ((X-Xm)**2 +(Y-Ym)**2 +(Z-Zm)**2)**1.5) ) +  (Ym/((Xm**2+Ym**2+Zm**2)**1.5)) )
fMoonZ = -GMM*( ( (Z-Zm)/ ( ((X-Xm)**2 +(Y-Ym)**2 +(Z-Zm)**2)**1.5) ) +  (Zm/((Xm**2+Ym**2+Zm**2)**1.5)) )

#Solar radiation pressure terms
fSRPX = Aratio*( (pSRP*(aS**2)*(X-Zs))/(((X-Xs)**2 +(Y-Ys)**2 +(Z-Zs)**2)**1.5) )
fSRPY = Aratio*( (pSRP*(aS**2)*(Y-Ys))/(((X-Xs)**2 +(Y-Ys)**2 +(Z-Zs)**2)**1.5) )
fSRPZ = Aratio*( (pSRP*(aS**2)*(Z-Zs))/(((X-Xs)**2 +(Y-Ys)**2 +(Z-Zs)**2)**1.5) )

#Combining final terms
d2X = fKepX# + fJ2X + fC22X + fS22X + fMoonX + fSunX + fSRPX
d2Y = fKepY# + fJ2Y + fC22Y + fS22Y + fMoonY + fSunY + fSRPY
d2Z = fKepZ# + fJ2Z + fC22Z + fS22Z + fMoonZ + fSunZ + fSRPZ

#Defining 1st order ODE system substitute variables
dX = pX
dY = pY
dZ = pZ

#Creating system of equations for integration
sysODE = [(X, dX), (Y, dY), (Z, dZ), (dX, d2X), (dY, d2Y), (dZ,d2Z)]

#Defining initial conditions
sysInit = [2.6533e7, 0, 0, 0, 2.2242e3, 3.1765e3]

#Creating integrator object
ta = hy.taylor_adaptive(sysODE,sysInit)

#Printing taylor object properties to stdout
print(ta)

#timestep grid
grid = np.linspace(0,31557600,100000)
#propagating along timesteps of grid
out = ta.propagate_grid(grid)
print(ta)

#saving propagator output
np.savetxt("propOut.txt", out[4],delimiter=',')
#saving time grid
np.savetxt("timeOut.txt", grid[:])
posMagnitude = np.sqrt((out[4][:, 3]**2) + (out[4][:, 4]**2) + (out[4][:, 5]**2))

import matplotlib.pyplot as plt

#3d positional orbit plot
fig = plt.figure(figsize=(9, 9))
ax = plt.axes(projection ='3d')

plt.plot(out[4][:, 3], out[4][:, 4],out[4][:, 5])

plt.grid()
plt.show()

#plot of axial position vs time
fig2 = plt.figure()
#plt.plot( grid[:], out[4][:, 3] - out[4][0, 3])
#plt.plot( grid[:], out[4][:, 4] - out[4][0, 4])
plt.plot( grid[:], posMagnitude[:])
plt.show()


#computing and printing script execution time
scriptExecTime = (time.time() - scriptStartTime)
print('Script Execution time: ' + str(scriptExecTime)[0:6] + ' seconds')