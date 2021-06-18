#This file contains coordinate system conversion functions
#For converting between a geocentric cartesian system and Orbital elements and opposite
from numpy.ma.core import dot

#importing NumPy library
import numpy as np

def vec2orb(pos, vel):
#Function for converting a state vector to orbital elements (Around Earth/Geocentric frame)
#UNITS - KM, KM/s (INPUTS)
#OUTPUT UNITS ARE: KM, RADIANS

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

#returning orbital elements
    return a, emag, i, raan, argp, nu


def orb2vec(a, e, i, raan, argp, nu):
#Function for converting from orbital elements to state vector
#INPUTS - semimajor axis-a, eccentricity-e, inclination-i, ascending node-raan
#       - argument of periapsis-argp, mean anomaly-M
#OUTPUTS
#State vector - 2 1x3 arrays of positon and velocity respectively
#UNITS KM, Radians, KM/s 

#importing NumPy library
    import numpy as np

#defining constants
    u = np.double(3.986004407799724e+5)   #GM Earth - km3/sec2

#calculating semilatus rectum
    sl = a*(1-e**2)

#calculating magnitude (polar coordinate system)
    rm = sl/(1+e*np.cos(nu))

#calculating angular position (polar coordinate system)
    lat = argp+nu

#calculating and combining position vector
    r = [(rm*(np.cos(raan)*np.cos(lat) - np.sin(raan)*np.cos(i)*np.sin(lat))), (rm*(np.sin(raan)*np.cos(lat) - np.cos(raan)*np.cos(i)*np.sin(lat))), (rm*np.sin(i)*np.sin(lat))]

#calculating velocity vector
    vx = -np.sqrt(u/sl) * (np.cos(raan)*(e* np.sin(argp) + np.sin(lat)) + np.sin(raan)*np.cos(i)*(e*np.cos(argp)+np.cos(lat)) )
    vy = -np.sqrt(u/sl) * (np.sin(raan)*(e* np.sin(argp) + np.sin(lat)) - np.cos(raan)*np.cos(i)*(e*np.cos(argp)+np.cos(lat)) )
    vz = np.sqrt(u/sl) * (e * np.cos(argp) + np.cos(lat)) * np.sin(i)
#combining velocity vector
    v = [vx, vy, vz]

#returning position and velocity vectors
    return r, v

#callback Function for integrator
#for use as a NON-TERMINAL event with the integrator

def callback1(ta, time, dsign):
    print("callback function triggered")



#Function for taking in the debris orbital elements

def ElementIn(elementFilePath,fileExtention):

    from os.path import exists

    limiter = " "
    elements = []
    i = 0

    while exists((elementFilePath + (str(i)) + fileExtention)):

        elementLine = np.loadtxt(fname = (elementFilePath + str(i) + fileExtention), delimiter = limiter)

        elements.append(list(elementLine))

        i += 1

    return(elements)
