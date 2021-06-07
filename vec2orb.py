from numpy.ma.core import dot
#Function for converting a state vector to orbital elements (Around Earth/Geocentric frame)
#UNITS - KM, KM/s (INPUTS)
#OUTPUT UNITS ARE: KM, RADIANS
def vec2orb(pos,vel):

#importing NumPy library
    import numpy as np

#Defining constants
    u = np.double(3.986004407799724e+5)   #GM Earth - km3/sec2

#calculating magnitude of position vector
    mag = np.linalg.norm(pos)
#calculating orbital momentum vector
    h = np.cross(pos,vel)

#calculating eccentricity vector
    e = (np.cross(vel,h)/u) - (pos/mag)

#calculating vector that points towards ascending node
    n = [-h[1], h[0], 0]

#calculating true anomaly (radians)
    nu = np.arccos(np.dot(e,pos)/(np.linalg.norm(e)*mag))

#true anomaly quadrant check
    if np.dot(pos,vel) < 0:
        nu = 2*np.pi - nu

#calculating inclination
    i = np.arccos(h[2]/np.linalg.norm(h))

#calculating eccentricity magnitude
    emag = np.linalg.norm(e)

#calculating eccentric anomaly
    E = 2*np.arctan(np.tan(nu/2)/np.sqrt((1+emag)/(1-emag)))

#calculating longitude of the ascending node
    raan = np.arccos(n[0]/np.linalg.norm(n))

#raan quadrant check
    if n[0] < 0:
        raan = 2*np.pi - raan

#calculating argument of periapsis
    argp = np.arccos(np.dot(n,e)/(np.linalg.norm(n)*emag))

#argp quadrant check
    if e[2] < 0:
        argp = 2*np.pi - argp
    
#calculating mean anomaly
    M = E - emag*np.sin(E)

#calculating semimajor axis
    a = 1/((2/mag)-((np.linalg.norm(vel)**2)/u))

    return a, emag, i, raan, argp, nu