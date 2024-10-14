##
## Algorithm from:
## Orbital Mechanics for Engineering Students, 
## 

import numpy as np
import math
from scipy.spatial.transform import Rotation as R

def elements2state(elements_dict, mu=398600.435507e9):
# get elements
    theta = elements_dict['true anomaly']
    h = elements_dict['specific angular momentum']
    e = elements_dict['eccentricity']
    Omega = elements_dict['Omega']
    omega = elements_dict['omega']
    i = elements_dict['inclination']
# trigonometric numbers
    st = np.sin(theta)
    ct = np.cos(theta)
# calculate position vector in perifocal frame
    r_pfcl = h**2/mu/(1e0+e*ct) * np.array((ct, st, 0e0))
# calculate velocity vector in perifocal frame
    v_pfcl = mu/h * np.array([-st, e+ct, 0e0])
# transformation matrix R, from perifocal to geocentric equatorial coordinates
    rot = R.from_euler('ZXZ', [Omega, i, omega])
# rotate and return the vectors, i.e. pos and vel, and scale if needed
    pos = rot.apply(r_pfcl)
    vel = rot.apply(v_pfcl)
    return pos, vel

def state2elements(pos, vel, mu=398600.435507e9):
# distance
    r = np.linalg.norm(pos)
# speed
    v = np.linalg.norm(vel)
# radial velocity
# Note that if vr > 0, the satellite is flying away from perigee. 
# If vr < 0, it is flying towards perigee.
    vr = np.dot(pos/r, vel)
# specific angular momentum and its magnitude [km^2/sec]
    h = np.cross(pos, vel)
    hn = np.linalg.norm(h)
# Inclination (i)
# Recall that i must lie between 0 and 180 [deg], so there is no quadrant 
# ambiguity. If 90 < i ≤ 180, the orbit is retrograde.
    inclination = np.arccos(h[2]/hn)
# vector defining the nodal line (and its magnitude)
    N = np.cross(np.array([0e0,0e0,1e0]), h)
    Nn = np.linalg.norm(N)
# Right Ascencion of the Ascending Node (Ω)
    Omega = 2 * np.pi - np.arccos(N[0] / Nn)
# eccentricity vector and magnitude
    e = np.cross(vel, h) / mu - pos / r
    eccentricity = np.linalg.norm(e)
# argument of perigee
    omega = 2 * np.pi - np.arccos(np.dot(N/Nn, e/eccentricity))
# true anomaly
    theta = np.arccos(np.dot(e/eccentricity,pos/r))
# a couple more quantities:
    perigee = (hn*hn/mu) * (1e0/(1+eccentricity))
    apogee  = (hn*hn/mu) * (1e0/(1-eccentricity))
    semimajor = .5*(perigee + apogee)
    period  = 2 * np.pi * semimajor * np.sqrt(semimajor) / (np.sqrt(mu)) # sec

    return {'specific angular momentum': hn,
            'inclination': inclination,
            'Omega': Omega,
            'eccentricity': eccentricity,
            'omega': omega,
            'true anomaly': theta,
            'perigee radius': perigee,
            'apogee radius': apogee,
            'semimajor': semimajor,
            'period': period}
