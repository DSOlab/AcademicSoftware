##
## Algorithm from:
## Orbital Mechanics for Engineering Students, 
## 

import numpy as np
import math
from scipy.spatial.transform import Rotation as R

def elements2state(elements_dict, mu=398600.435507e9, useKm=False):
# transform all meters to km
    mu = mu * 1e-9 ## [km^3*s^2]
# get elements
    theta = elements_dict['true anomaly']
    h = elements_dict['specific angular momentum']
    e = elements_dict['eccentricity']
    Omega = elements_dict['Omega']
    omega = elements_dict['omega']
    i = elements_dict['inclination']
# trigonometric numbers
    st = math.sin(theta)
    ct = math.cos(theta)
# calculate position vector in perifocal frame
    r_pfcl = (h*h/mu) * (1e0/(1e0+e*ct)) * np.array([ct, st, 0e0])
# calculate velocity vector in perifocal frame
    v_pfcl = (mu/h) * np.array([-st, e+ct, 0e0])
# transformation matrix R, from perifocal to geocentric equatorial coordinates
    rot= R.from_euler('zxz', [omega, i, Omega])
# rotate and return the vectors, i.e. pos and vel, and scale if needed
    scale = 1e0 if useKm else 1e3
    return rot.apply(r_pfcl) * scale, rot.apply(v_pfcl) * scale

def state2elements(pos, vel, mu=398600.435507e9, useKm=False):
# transform all meters to km
    mu = mu * 1e-9 ## [km^3*s^2]
    r = pos * 1e-3 ## km
    v = vel * 1e-3 ## km/sec
# distance
    rn = np.linalg.norm(r)
# speed
    vn = np.linalg.norm(v)
# radial velocity
# Note that if vr > 0, the satellite is flying away from perigee. If vr < 0, 
# it is flying towards perigee.
    vr = np.dot(r,v) / rn
# angular momentum and its magnitude
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
# Inclination (i)
# Recall that i must lie between 0 and 180 [deg], so there is no quadrant 
# ambiguity. If 90 < i ≤ 180, the orbit is retrograde.
    inclination = math.acos(h[2]/hn)
# vector defining the nodal line (and its magnitude)
    N = np.cross(np.array([0,0,1]), h)
    Nn = np.linalg.norm(N)
# Right Ascencion of the Ascending Node (Ω)
    Omega = math.acos(N[0]/Nn) if N[1]>=0 else math.pi*2 - math.acos(N[0]/Nn)
# eccentricity vector and magnitude
    e = (1e0/mu) * (np.cross(v,h) - (mu/rn)*r)
    eccentricity = np.linalg.norm(e)
# argument of perigee
    omega = math.acos(np.dot(N,e) / (Nn*eccentricity));
    if e[2] < 0: omega = math.pi * 2 - omega
# true anomaly
    theta = np.dot(e,r) / (eccentricity*rn)
    if vr < 0: theta = math.pi * 2 - theta

# acouple more quantities:
    perigee = (hn*hn/mu) * (1e0/(1+eccentricity))
    apogee  = (hn*hn/mu) * (1e0/(1-eccentricity))
    semimajor = .5*(perigee + apogee)
    period  = 2 * math.pi * semimajor * math.sqrt(semimajor) / (math.sqrt(mu))

    # check units before returning
    if not useKm: 
        hn *= 1e6
        perigee *= 1e3
        apogee *= 1e3
        semimajor *= 1e3

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
